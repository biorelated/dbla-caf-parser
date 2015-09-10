require 'bio'
require File.join(File.expand_path(File.dirname(__FILE__)), 'dbla-caf-parser','caf/caf_parser')
#require 'progress_bar'

caf_file = "#{ENV['HOME']}/Batch1/batch1_cafs/batch1.caf"
parser = Parser.new(caf_file)

parser.reads.each do |read|
 #puts read.to_fasta(:raw)
 #puts read.to_fasta(:without_vector)
 puts read.qual_str
end
