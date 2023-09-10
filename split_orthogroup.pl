#!/usr/bin/perl -w

# Find all lines in which the pattern (\w+\|) is matched one or more times.
# If all matches are of the same string, write to file 2.  If any matches are
# not the same, write to file 1.

# Comment this out in order to put the filename on the command line.
@ARGV = ("orthogroups_edited.txt");

# Change to the filenames you want to use.  Or write one to STDOUT and the
# other to STDERR.
open(FILE1, "> real_orthogroups_perl.txt") or die "Can't open 'real_orthogroups_perl.txt'\n$!";
open(FILE2, "> singleton_and_inparalogs_perl.txt") or die "Can't open 'singleton_and_inparalogs_perl.txt'\n$!";

while (<>) {
    if (m# \w+\|#g) {
        $found = $&;
        # print "Starting search for '$found'\n";
        $total = 0;
        $same = 0;
        while (m# \w+\|#g) {
            $total++;
            $same++ if $found eq $&;
        }
        # print "total ($total) == same ($same)?\n";
        if ($total == $same) {
            print FILE2 $_;
        } else {
            print FILE1 $_;
        }
    }
}
close(FILE1);
close(FILE2);
exit 0;