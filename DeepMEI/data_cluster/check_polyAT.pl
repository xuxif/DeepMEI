#!/usr/bin/perl
use strict;
use warnings;

# Function to check for poly-A/T sequence in a read
sub check_polyAT {
    my ($cigar, $seq) = @_;
    my ($soft_clip_seq, $match_seq);

    # Check for soft-clipping at the start of the sequence
    if ($cigar =~ /^(\d+)S/) {
        my $soft_clip_length = $1;
        $soft_clip_seq = reverse(substr($seq, 0, $soft_clip_length));
        $match_seq = substr($seq, $soft_clip_length,15);
    }
    # Check for soft-clipping at the end of the sequence
    elsif ($cigar =~ /(\d+)S$/) {
        my $soft_clip_length = $1;
        $soft_clip_seq = substr($seq, -$soft_clip_length);
        $match_seq = reverse(substr($seq,-1*($soft_clip_length+15) , 15));
    } else {
        return 0;  # No soft-clipping found
    }

    # Check for poly-A or poly-T sequence in the soft-clipped part
    #    if ($soft_clip_seq =~ /^A{3,}/ or $soft_clip_seq=~/^T{3,}/) {
        # Check for a longer poly-A or poly-T in the matched part
        if ($match_seq =~ /^[A]{10,}/ or $match_seq =~/^[T]{10,}/) {
            return 1;  # Match found
	    #    }
    }
    return 0;  # No match found
}

# Process the SAM file
sub process_file {
    my ($filename, $record) = @_;
    open my $fp, "<", $filename or die "Cannot open file $filename: $!";

    while (my $line = <$fp>) {
        next if $line =~ /^\@/;  # Skip header lines
        chomp $line;
        my @fields = split(/\t/, $line);
        my $seq = $fields[9];
        my $cigar = $fields[5];

        if (check_polyAT($cigar, $seq)) {
            my $formatted_record = $record;
            $formatted_record =~ s/[:\-]/\t/g;
            return 1;
            last;
        }
    }
    close $fp;
    return 0;
}

my $record = $ARGV[0];
my $left=process_file("HG002_${record}_mapClipL.sam", $record);
my $right=process_file("HG002_${record}_mapClipR.sam", $record);
if($left*$right==1)
{
	$record=~s/[:\-]/\t/g;
	print "$record\n";
}
