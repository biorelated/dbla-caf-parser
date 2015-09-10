require File.join(File.expand_path(File.dirname(__FILE__)),'read')

class Parser

  def initialize(file)
    @file = file
  end

  def reads
    $/ = "\n\n"
    records = []
    File.open(@file) do |f|
      f.each_line do |line|
        read = Read.new
        case line
        when /^DNA/
          records << Read.new
          records.last.dna_data = line.strip
        when /^BaseQuality/
          records.last.quality_data = line.strip
        when /^Sequence/
          records.last.other_data = line.strip
        else
          puts "Unrecognized line: #{line}"
        end
      end
    end
    records
  end
end
