#!/usr/bin/perl -w
use strict;

my $line = "";
my @relations;
my $relationsFile = $ARGV[0];
my $dependenciesFile = $ARGV[1];
my $relationSetsFile = $ARGV[2];
my $processDependenciesExe = $ARGV[3];
my $calcroot = $ARGV[4];
my $outfileRoot = $relationsFile;
$outfileRoot =~ s/[\.].*$//;

sub edit_root_config_file($)
{
   my $outfile = shift;
   open ROOT_CONFIG, "<root.cfg" or die "Cannot open root.cfg : $!";
   open WORKFILE, ">root.cfg.work" or die "Cannot open root.cfg.work : $!";
   my $line = "";
   while ($line = <ROOT_CONFIG>)
   {
      chomp($line);
      if ($line =~ m/^RELATION_FILE/)
      {
         print WORKFILE "RELATION_FILE = $outfile\n";
      }
      else
      {
         print WORKFILE "$line\n";
      }
   }
   close ROOT_CONFIG;
   close WORKFILE;
   rename "root.cfg.work","root.cfg";
}

sub calc_root()
{
   my $rc = system("$calcroot");
   return $rc;
}

my $relationCount = `wc -l $relationsFile`;
$relationCount =~ s/ .*$//;
print STDERR "$relationsFile has $relationCount relations\n";
my $dependencyCount = `wc -l $dependenciesFile`;
$dependencyCount =~ s/ .*$//;
print STDERR "$dependenciesFile has $dependencyCount dependencies\n";

my $dep = 0;
for ($dep = 0; $dep < $dependencyCount; $dep++)
{
   my $cmd = "$processDependenciesExe -r $relationsFile -d $dependenciesFile -rsi $relationSetsFile -f $dep";
   print STDERR "Running [$cmd]\n";
   system($cmd);
   my $outfile = "$outfileRoot.rel$dep";
   edit_root_config_file($outfile);
   last if (calc_root() == 0);
}

exit(0);

my $rel = 0;
for ($rel = 0; $rel < $relationCount; $rel++)
{
   $relations[$rel]{ CHOSEN } = 0;
}

my @relationSets = ();
my @relationSetsCount = ();
open RELSETS, "<$relationSetsFile" or die "Cannot open $relationSetsFile : $!\n";
$line = <RELSETS>;
chomp($line);
my $relationSetCount = int($line);
my $i = 0;
while ($line = <RELSETS>)
{
   chomp($line);
   my ($relcount, @rels) = split(' ', $line);
   my $rel;
   my $count = 0;
   foreach $rel (@rels)
   {
      $relationSets[$i][$count] = int($rel);
      $count++;
   }
   $relationSetsCount[$i] = $count;
   if ($relcount != $count)
   {
      print STDERR "Problem: relcount = $relcount, count = $count\n";
   }
   $i++;
}
close RELSETS;
print STDERR "$i relation sets read from $relationSetsFile\n";

my $min_ones = 1000000000;
my @best_deps = ();
my $best_line = "";
open DEPENDENCIES, "<$dependenciesFile" or die "Cannot open $dependenciesFile : $!\n";
my $j = 0;
my $best_j = 0;
while ($line = <DEPENDENCIES>)
{
   chomp($line);
   my @deps = split('',$line);
   my $ones = 0;
   my $dep = 0;
   my $i = 0;
   my $chosen = 0;
   foreach $dep (@deps)
   {
      if ("$dep" eq "1")
      {
         $ones++;
	 # choose relation set $i
	 my $j;
	 for ($j = 0; $j < $relationSetsCount[$i]; $j++)
	 {
            my $rel = $relationSets[$i][$j];
	    if ($relations[$rel]{ CHOSEN } == 0)
	    {
	       $relations[$rel]{ CHOSEN } = 1;
	       $chosen++;
	    }
	    else
	    {
	       $relations[$rel]{ CHOSEN } = 0;
	       $chosen--;
	    }
         }
      }
      $i++;
   }
   if ($chosen == 0 or $chosen % 2 != 0)
   {
      for ($i = 0; $i < $relationCount; $i++)
      {
         $relations[$i]{ CHOSEN } = 0;
      }
      next;
   }
   my $outfile = "$outfileRoot.rel$j";
   my $outfile1 = "$outfileRoot.rel$j.rat";
   open OUTFILE, ">$outfile" or die "Cannot open $outfile : $!\n";
   open OUTFILE1, ">$outfile1" or die "Cannot open $outfile1 : $!\n";
   open RELATIONS, "<$relationsFile" or die "Cannot open $relationsFile : $!\n";

   for ($i = 0; $i < $relationCount; $i++)
   {
      my $line = <RELATIONS>;
      if ($relations[$i]{ CHOSEN } == 1)
      {
	 chomp($line);
         my ($a, $b, $fac1, $fac2) = ($line =~ m/^([^ ]*) ([^ ]*) : ([^:]*) : ([^:]*) :.*$/);
         $fac1 =~ s/\/[0-9]+ / /g;
         $fac1 =~ s/\/[0-9]+$//g;
         print OUTFILE "$a $b $fac1\n";
         print OUTFILE1 "$a $b $fac2\n";
         $relations[$i]{ CHOSEN } = 0;
      }
   }
   print OUTFILE "!\n";
   print OUTFILE1 "!\n";
   close OUTFILE;
   close OUTFILE1;
   close RELATIONS;
   if ($chosen % 2 == 0)
   {
      if ($chosen > 0 and $chosen < $min_ones)
      {
         $min_ones = $chosen;
         @best_deps = @deps;
         $best_line = $line;
         $best_j = $j;
      }
      $j++;
   }
   edit_root_config_file($outfile);
   last if (calc_root() == 0);
}
close DEPENDENCIES;

print "$j dependencies found\n";
print "Best dependency ($best_j) which has $min_ones relations\n";

