#!/usr/bin/perl -w
use strict;

my $bindir = "c:/maths/gnfs/vcbin";
my $filter = $bindir . "/filter";
my $buildMatrix = $bindir . "/buildMatrix";
my $blockLanczos = $bindir . "/bl";
my $processDependenciesPl = $bindir . "/processDependencies.pl";
my $processDependenciesExe = $bindir . "/processDependencies";
my $calcroot = $bindir . "/calcroot";
my $sort = "/bin/sort";

my $workfile1 = "workfile1";
my $workfile2 = "workfile2";

# variables set from command line
my $input_relation_file;
my $output_relation_file;
my $relation_set_file;
my $output_primes_file = 0;
my $relation_set_primes_file;
my $relation_set_infile;
my $input_primes_file = 0;
my $relation_set_primes_infile;
my $matrix;
my $checkpointfile;
my $dependencies;
my $merge_level;
my $min_merge_level;
my $merge_only = 0;
my $max_passes;
my $filtmin;
my $excess;

my $command = 0; 
my $command_count = 0;

my $REMOVE_DUPLICATES = 1;
my $REMOVE_SINGLETONS = 2;
my $PRUNE = 4;
my $MERGE = 8;
my $BUILD_MATRIX = 16;
my $SOLVE_MATRIX = 32;
my $PROCESS_DEPENDENCIES = 64;
my $CALC_ROOT = 128;

my %commands = ( 
   "REMOVE_DUPLICATES" => $REMOVE_DUPLICATES,
   "R" => $REMOVE_DUPLICATES,
   "REMOVE_SINGLETONS" => $REMOVE_SINGLETONS,
   "RS" => $REMOVE_SINGLETONS,
   "PRUNE" => $PRUNE,
   "P" => $PRUNE,
   "MERGE" => $MERGE,
   "M" => $MERGE,
   "BUILD_MATRIX" => $BUILD_MATRIX,
   "B" => $BUILD_MATRIX,
   "SOLVE_MATRIX" => $SOLVE_MATRIX,
   "S" => $SOLVE_MATRIX,
   "CALC_ROOT" => $CALC_ROOT,
   "C" => $CALC_ROOT,
   "PROCESS_DEPENDENCIES" => $PROCESS_DEPENDENCIES,
   "D" => $PROCESS_DEPENDENCIES,
   "MERGE-BUILD_MATRIX" => ($MERGE + $BUILD_MATRIX),
   "M-B" => ($MERGE + $BUILD_MATRIX),
   "MERGE-SOLVE_MATRIX" => ($MERGE + $BUILD_MATRIX + $SOLVE_MATRIX),
   "M-S" => ($MERGE + $BUILD_MATRIX + $SOLVE_MATRIX),
   "MERGE-PROCESS_DEPENDENCIES" => ($MERGE + $BUILD_MATRIX + $SOLVE_MATRIX + $PROCESS_DEPENDENCIES),
   "M-D" => ($MERGE + $BUILD_MATRIX + $SOLVE_MATRIX + $PROCESS_DEPENDENCIES),
   "BUILD_MATRIX-SOLVE_MATRIX" => ($BUILD_MATRIX + $SOLVE_MATRIX),
   "B-S" => ($BUILD_MATRIX + $SOLVE_MATRIX),
   "BUILD_MATRIX-PROCESS_DEPENDENCIES" => ($BUILD_MATRIX + $SOLVE_MATRIX + $PROCESS_DEPENDENCIES),
   "B-D" => ($BUILD_MATRIX + $SOLVE_MATRIX + $PROCESS_DEPENDENCIES),
   "SOLVE_MATRIX-PROCESS_DEPENDENCIES" => ($SOLVE_MATRIX + $PROCESS_DEPENDENCIES), 
   "S-D" => ($SOLVE_MATRIX + $PROCESS_DEPENDENCIES), 
   );

sub command_includes($)
{
   my $c = shift;
   my $rc = ($command & $c);
   return $rc;
}

sub usage()
{
   print STDERR "Usage : postProcessing.pl -r relations\n";
   print STDERR "                         [-rs relation_sets]\n";
   print STDERR "                         [-rsp [relation_sets_primes] ]\n";
   print STDERR "                         [-ro relations_output]\n";
   print STDERR "                         [-m matrix]\n";
   print STDERR "                         [-cp checkpointfile]\n";
   print STDERR "                         [-d dependencies]\n";
   print STDERR "                         [-ml merge_level]\n";
   print STDERR "                         [-mm min_merge_level]\n";
   print STDERR "                         [-merge_only]\n";
   print STDERR "                         [-rsi input_relation_sets]\n";
   print STDERR "                         [-rspi [input_relation_sets_primes] ]\n";
   print STDERR "                         [-mp max_passes]\n";
   print STDERR "                         [-f filtmin]\n";
   print STDERR "                         [-e excess]\n";
   print STDERR "                         -c command[-command] ]\n";
   print STDERR "                         [-h[elp]]\n";
   print STDERR "where command is one of \n";
   print STDERR "   R[EMOVE_DUPLICATES]\n";
   print STDERR "   [RS|REMOVE_SINGLETONS]\n";
   print STDERR "   P[RUNE]\n";
   print STDERR "   M[ERGE]\n";
   print STDERR "   B[UILD_MATRIX]\n";
   print STDERR "   S[OLVE_MATRIX]\n";
   print STDERR "   [D|PROCESS_DEPENDENCIES]\n";
   print STDERR "   C[ALC_ROOT]\n";
   exit(1);
}

sub err_usage($)
{
   my $message = shift;
   print STDERR "$message\n";
   usage();
}

sub check_command_line()
{
   if (command_includes($REMOVE_DUPLICATES))
   {
      err_usage("input relation file not supplied with -r option") if (not defined $input_relation_file);
      err_usage("output relation file not supplied with -ro option") if (not defined $output_relation_file);
      $merge_level = 0 if (not defined $merge_level);
      return;
   }

   if (command_includes($REMOVE_SINGLETONS))
   {
      err_usage("input relation file not supplied with -r option") if (not defined $input_relation_file);
      err_usage("output relation file not supplied with -ro option") if (not defined $output_relation_file);
      $merge_level = 1;
      $filtmin = 1000 if (not defined $filtmin);
      $excess = 500 if (not defined $excess);
      return;
   }

   if (command_includes($PRUNE))
   {
      err_usage("input relation file not supplied with -r option") if (not defined $input_relation_file);
      err_usage("output relation file not supplied with -ro option") if (not defined $output_relation_file);
      $merge_level = 1 if (not defined $merge_level);
      $filtmin = 1000 if (not defined $filtmin);
      $excess = 500 if (not defined $excess);
      return;
   }

   if (command_includes($MERGE))
   {
      $merge_level = 2 if (not defined $merge_level);
      $filtmin = 1000 if (not defined $filtmin);
      $excess = 500 if (not defined $excess);
      # input files:
      #    input_relation_file 
      #    relation_set_infile
      #    relation_set_primes_infile
      #
      $relation_set_primes_infile = $relation_set_infile . ".primes" if (not defined $relation_set_primes_infile and $input_primes_file);
      if ($merge_only)
      {
         #err_usage("must supply input relation sets file with -rsi") if (not defined $relation_set_infile);
      }
      else
      {
         err_usage("must supply input relation file with -r") if (not defined $input_relation_file);
      }

      # possible output files:
      #    output_relation_file
      #    relation_set_file
      #    relation_set_primes_file
      $relation_set_file = $output_relation_file . ".sets" if (not defined $relation_set_file and defined $output_relation_file);
      $relation_set_primes_file = $relation_set_file . ".primes" if (not defined $relation_set_primes_file and defined $relation_set_file and $output_primes_file);
      if ($merge_only)
      {
         err_usage("must supply output relation sets file with -rs") if (not defined $relation_set_file);
      }
      else
      {
         err_usage("must supply output relation file with -ro") if (not defined $output_relation_file);
      }
      $max_passes = 2 if (not defined $max_passes);
   }

   if (command_includes($BUILD_MATRIX))
   {
      err_usage("input relation file not supplied with -r option") if (not defined $input_relation_file);
      err_usage("matrix file name not supplied with -m option") if (not defined $matrix);
      err_usage("filtmin not supplied with -f option") if (not defined $filtmin);
      $relation_set_infile = $input_relation_file . ".sets" if (not defined $relation_set_infile);
      $relation_set_primes_infile = $relation_set_infile . ".primes" if (not defined $relation_set_primes_infile and defined $relation_set_infile);
   }

   if (command_includes($SOLVE_MATRIX))
   {
      err_usage("matrix file name not supplied with -m option") if (not defined $matrix);
      err_usage("dependencies file name not supplied with -d option") if (not defined $dependencies);
   }

   if (command_includes($PROCESS_DEPENDENCIES))
   {
      err_usage("input relation file not supplied with -r option") if (not defined $input_relation_file);
      err_usage("dependencies file name not supplied with -d option") if (not defined $dependencies);
      $relation_set_infile = $input_relation_file . ".sets" if (not defined $relation_set_infile);
      return;
   }

   if (command_includes($CALC_ROOT))
   {
      return;
   }
}

sub process_command_line()
{
   my $argc = @ARGV;
   my $arg = 0;
   while ($arg < $argc)
   {
      if ($ARGV[$arg] eq "-r")
      {
         $arg++;
	 $input_relation_file = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-rs")
      {
         $arg++;
	 $relation_set_file = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-rsi")
      {
         $arg++;
	 $relation_set_infile = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-rsp")
      {
         $output_primes_file = 1;
         if ($ARGV[$arg + 1] !~ "^-")
	 {
            $arg++;
	    $relation_set_primes_file = $ARGV[$arg];
         }
      }
      elsif ($ARGV[$arg] eq "-rspi")
      {
         $input_primes_file = 1;
         if ($ARGV[$arg + 1] !~ "^-")
	 {
            $arg++;
	    $relation_set_primes_infile = $ARGV[$arg];
         }
      }
      elsif ($ARGV[$arg] eq "-ro")
      {
         $arg++;
	 $output_relation_file = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-m")
      {
         $arg++;
	 $matrix = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-cp")
      {
	 $arg++;
	 $checkpointfile = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-d")
      {
         $arg++;
	 $dependencies = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-ml")
      {
         $arg++;
	 $merge_level = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-mm")
      {
         $arg++;
	 $min_merge_level = $ARGV[$arg];
	 $min_merge_level = 3 if ($min_merge_level < 3);
      }
      elsif ($ARGV[$arg] eq "-merge_only")
      {
         $merge_only = 1;
      }
      elsif ($ARGV[$arg] eq "-mp")
      {
         $arg++;
	 $max_passes = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-f")
      {
         $arg++;
	 $filtmin = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-e")
      {
         $arg++;
	 $excess = $ARGV[$arg];
      }
      elsif ($ARGV[$arg] eq "-c")
      {
         $arg++;
	 my $c = $ARGV[$arg];
	 usage() if (not defined $c);
	 usage() if (not exists $commands { $c });
	 $command = $commands { $c };
      }
      elsif ($ARGV[$arg] eq "-h" or $ARGV[$arg] eq "-help")
      {
         usage();
      }
      $arg++;
   }
   check_command_line();
}

sub remove_duplicates($$)
{
   my $input_relation_file = shift;
   my $output_relation_file = shift;
   print "$filter -r $input_relation_file -m 0 -ro $output_relation_file\n";
   system("$filter -r $input_relation_file -m 0 -ro $output_relation_file");
}

sub remove_singletons($$$$)
{
   my $input_relation_file = shift;
   my $output_relation_file = shift;
   my $filtmin = shift;
   my $excess = shift;
   print "$filter -r $input_relation_file -m 1 -f $filtmin -e $excess -ro $output_relation_file\n";
   system("$filter -r $input_relation_file -m 1 -f $filtmin -e $excess -ro $output_relation_file");
}

sub prune($$$$)
{
   my $input_relation_file = shift;
   my $output_relation_file = shift;
   my $filtmin = shift;
   my $excess = shift;
   print "$filter -r $input_relation_file -m 2 -f $filtmin -e $excess -ro $output_relation_file\n";
   system("$filter -r $input_relation_file -m 2 -f $filtmin -e $excess -ro $output_relation_file");
}

sub merge($$$$$$$)
{
   my $input_relation_file = shift;
   my $output_relation_file = shift;
   my $relation_set_file = shift;
   my $merge_level = shift;
   my $filtmin = shift;
   my $excess = shift;
   my $max_passes = shift;
   my $cmd = "$filter";
   $cmd .= " -r $input_relation_file" if (defined $input_relation_file);
   $cmd .= " -ro $output_relation_file" if (defined $output_relation_file);
   $cmd .= " -rspo $relation_set_primes_file" if (defined $relation_set_primes_file);
   if ($merge_only)
   {
      $cmd .= " -merge_only";
      $cmd .= " -rsi $relation_set_infile" if (defined $relation_set_infile);
      $cmd .= " -rspi $relation_set_primes_infile" if (defined $relation_set_primes_infile);
   }
   $cmd .= " -m $merge_level";
   $cmd .= " -mm $min_merge_level" if (defined $min_merge_level);
   $cmd .= " -f $filtmin -e $excess -mp $max_passes -rso $relation_set_file";
   print "[$cmd]\n";
   system("$cmd");
}

sub build_matrix($$$$$)
{
   my $input_relation_file = shift;
   my $relation_set_file = shift;
   my $relation_set_prime_file = shift;
   my $matrix = shift;
   my $filtmin = shift;
   print "$buildMatrix -f $filtmin -r $input_relation_file -rs $relation_set_file -rsp $relation_set_prime_file -m $matrix -c sieve.cfg\n";
   system("$buildMatrix -f $filtmin -r $input_relation_file -rs $relation_set_file -rsp $relation_set_prime_file -m $matrix -c sieve.cfg");
}

sub solve_matrix($$$)
{
   my $matrix = shift;
   my $dependencies = shift;
   my $checkpointfile = shift;
   my $workfile = $dependencies . ".work";
   my $cmd = "$blockLanczos -split -m $matrix";
   $cmd .= " -c $checkpointfile" if (defined $checkpointfile);
   $cmd .= " > $workfile";
   print "[$cmd]\n";
   system("$cmd");
   system("$sort -u $workfile > $dependencies");
   #unlink $workfile;
}

sub process_dependencies($$$)
{
   my $input_relation_file = shift;
   my $relation_set_file = shift;
   my $dependencies = shift;
   system("$processDependenciesPl $input_relation_file $dependencies $relation_set_file $processDependenciesExe $calcroot");
}

sub calc_root()
{
   system("$calcroot");
}

sub print_name_value_pair($$)
{
   my $name = shift;
   my $value = shift;
   print STDERR "$name [";
   print STDERR "$value" if (defined $value);
   print STDERR "]\n";
}

sub process_commands()
{
   if (command_includes($REMOVE_DUPLICATES))
   {
      print STDERR "Removing duplicates relations from $input_relation_file, output to $output_relation_file\n";
      remove_duplicates($input_relation_file, $output_relation_file);
   }
   elsif (command_includes($REMOVE_SINGLETONS))
   {
      print STDERR "Removing singleton relation from $input_relation_file, output to $output_relation_file\n";
      remove_singletons($input_relation_file, $output_relation_file, $filtmin, $excess);
   }
   elsif (command_includes($PRUNE))
   {
      print STDERR "Pruning $input_relation_file, output to $output_relation_file\n";
      prune($input_relation_file, $output_relation_file, $filtmin, $excess);
   }
   elsif (command_includes($MERGE))
   {
      print STDERR "Merging:\n";
      print_name_value_pair("Input relation file", $input_relation_file);
      print_name_value_pair("Output relation file", $output_relation_file);
      print_name_value_pair("Input relation set file", $relation_set_infile);
      print_name_value_pair("Output relation file", $relation_set_file);
      print_name_value_pair("Input relation set primes file", $relation_set_primes_infile);
      print_name_value_pair("Output relation set primes file", $relation_set_primes_file);
      print_name_value_pair("Merge level", $merge_level);
      print_name_value_pair("Minimum merge level", $min_merge_level);
      print_name_value_pair("Minimum prime to filter", $filtmin);
      print_name_value_pair("Excess", $excess);
      print_name_value_pair("Maximum passes", $max_passes);
      merge($input_relation_file, $output_relation_file, $relation_set_file, $merge_level, $filtmin, $excess, $max_passes);
      $input_relation_file = $output_relation_file;
   }

   if (command_includes($BUILD_MATRIX))
   {
      print STDERR "Building matrix from relations in $input_relation_file, relation sets in $relation_set_infile, relation set primes in $relation_set_primes_infile and output to $matrix\n";
      build_matrix($input_relation_file, $relation_set_infile, $relation_set_primes_infile, $matrix, $filtmin);
   }

   if (command_includes($SOLVE_MATRIX))
   {
      print STDERR "Solving matrix in $matrix, output dependencies to $dependencies\n";
      solve_matrix($matrix, $dependencies, $checkpointfile);
   }

   if (command_includes($PROCESS_DEPENDENCIES))
   {
      print STDERR "Processing dependencies from $dependencies, relations in $input_relation_file, relation sets in $relation_set_infile\n";
      process_dependencies($input_relation_file, $relation_set_infile, $dependencies);
   }

   if (command_includes($CALC_ROOT))
   {
      print STDERR "Calculating square root\n";
      calc_root();
   }
}

process_command_line();
process_commands();

exit(0);
