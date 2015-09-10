class Read
  attr_accessor :dna_data,:quality_data,:other_data

  #extract the name for this read
  def name
    regexp = to_regexp("(VAR[PB|BP].*)$")
    dna_data.to_s.scan(regexp).join
  end

  def raw_seq
    regexp = to_regexp("(^DNA : VARPB.*)$")
    dna_data.to_s.gsub(regexp,'').delete("\n").strip
  end

  def raw_quals
    regexp = to_regexp("(^\\d.*)")
    quality_data.to_s.scan(regexp).to_s
  end

  def vector_pos
    regexp = to_regexp("^Seq_vec\\s+SVEC\\s+(\\d+\\s+\\d+)")
    find_vector_pos(regexp)
  end

  def clip_pos
    clipping = to_regexp("^Clipping\\s+QUAL\\s+(\\d+\\s+\\d+)")
    find_clip_pos(clipping)
  end

  def without_vector
    highest_start = vector_clips[0]
    raw_seq[highest_start, clip_len].to_s
  end

  # returns base qualities devoid of vector quality scores.
  def quals
    highest_start = vector_clips[0]
    raw_quals.split(/ /).slice(highest_start, clip_len).map(&:to_i) #if without_vector.size > 2 
  end

  # print base qualities 
  def qual_str
    ">#{name}\n #{quals.join(' ')}" if without_vector.size > 2
  end

  def vector_clips
    #first step: set highest_start to svec1_stop and lowest_end to 10 000
    highest_start = vector_pos[1]
    lowest_end = 10_000

    #second step
    #check whether svec2_start is less than 10 000; if true set lowest_end to svect2_start
    unless vector_pos[2] == 0
      lowest_end = vector_pos[2]  if vector_pos[2] < lowest_end
    else
      lowest_end = clip_pos[1]  #if we only have 1 svec entry,
    end
    highest_start = clip_pos[0] if highest_start < clip_pos[0]
    lowest_end    = clip_pos[1] if lowest_end > clip_pos[1]
    [highest_start, lowest_end]
  end

  #TODO this method is not accurate and is buggy..
  def without_primer
    begin
      seq = without_vector
      five_prime_primers.each do |primer|
        seq.gsub!(primer) {|x| 'X' * x.length}
      end

      #Replace all 5' bases before "X" with "X".
      seq.sub!(/\A[^X]+X/) { |x| "X" * x.length }

      three_prime_primers.each do |primer|
        seq.gsub!(primer) {|x| 'X' * x.length}
      end

      #Replace all 3' bases after "X" with "X".
      seq.sub!(/X[^X]+\z/) { |x| "X" * x.length }

      seq.gsub!(/X/,'')

      #start = seq.index(/[X]/).to_i
      #stop  = seq.rindex(/[X]/).to_i
      #clip = stop - start

      #unless start.nil? || stop.nil? || clip > 10
      #seq[start,clip + 1].gsub(/X/,'')
      #else
      #without_vector
      #end
    rescue
    end
  end

  def to_fasta(output=:raw)
    case output 
    when :without_primer
      Bio::Sequence.new(without_primer).output(:fasta,:header=>name,:width=>1000) if without_primer.size > 10
    when :without_vector
      Bio::Sequence.new(without_vector).output(:fasta,:header=>name,:width=>1000) if without_vector.size > 2
    when :raw
      Bio::Sequence.new(raw_seq).output(:fasta,:header=>name,:width=>1000)
    else
      puts "specifiy output format"
    end
  end

  private

  def to_regexp(string,opts={})
    options = {:case_sensitive => true}.merge!(opts)
    Regexp.new(string,options)
  end

  def primers
    ['G*CACG[A|C]AGTTT[C|T]GC','G*CCCATTC[G|C]TCGAACCA','GC[G|A]AAACT[T|G]CGTGC','TGGTTCGA[C|G]GAATGGGC'].map! { |primer| to_regexp(primer) }
  end

  def five_prime_primers
    ['G*CACG[A|C]AGTTT[C|T]GC','G*CCCATTC[G|C]TCGAACCA'].map! { |primer| to_regexp(primer) }
  end

  def three_prime_primers
    ['GC[G|A]AAACT[T|G]CGTGC','TGGTTCGA[C|G]GAATGGGC'].map! { |primer| to_regexp(primer) }
  end

  def find_vector_pos(regexp)
    positions                 = other_data.to_s.scan(regexp)
    svec1_startq, svec1_stopq = positions[0].to_s.split(/ /)
    svec2_startq, svec2_stopq = positions[1].to_s.split(/ /)
    [svec1_startq.to_i, svec1_stopq.to_i, svec2_startq.to_i, svec2_stopq.to_i]
  end

  def find_clip_pos(clipping)
    clipping_positions = other_data.to_s.scan(clipping)
    clip_qual_startq, clip_qual_stopq = clipping_positions[0].to_s.split(/ /)
    [clip_qual_startq.to_i, clip_qual_stopq.to_i]
  end

  def clip_len
    highest_start = vector_clips[0]
    lowest_end    = vector_clips[1]
    lowest_end.to_i - highest_start.to_i
  end

end
