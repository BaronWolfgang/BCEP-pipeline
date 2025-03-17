#!/usr/bin/perl
my $i=0;
my $command=$0;
while (defined($ARGV[$i])){
  $command .= " $ARGV[$i]";
  $i++;
}

use Getopt::Std;
use Cwd;
use lib "$ENV{NetSurfP}/modules";
use HOW;
use FASTA;
use BLAST;
use PROF;
#use strict;
#
# Default parameters
#
my $TYPE="HOW";
my $cmd;
#
# External programs
#
my $blastpgp = "$ENV{BLASTPROG}";

#
# Blast database
#
my $db = "nr";
my $itr = 4;
my $evalue = 0.00001;
my $tmpdir = "HOW$$";
my $alignments=0;
my $chunks=1;
#
# Process command line
#
getopts('hi:I:d:pP:S:lvkj:e:L:t:o:b:BRa:')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i name] [-I input_file_type] [-d name] [-j integer] [-e number] [-t name] [-k] [-L name] [-v] [-o name] \n");
  print ("Description:\n");
  print ("$0 - Make profiles for the how2004 neural network\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input file\n");
  print ("  -I  : input file type ([$TYPE]|FASTA|PROF)\n");
  print ("  -d  : Blast database [$db]\n");
  print ("  -j  : Number of psi-blast iterations [$itr]\n");
  print ("  -e  : psi-blast e-value [$evalue]\n");
  print ("  -t  : Temporary directory [$tmpdir]\n");
  print ("  -k  : Remove temporary directory [on]\n");
  print ("  -L  : Logfile [STDERR]\n");
  print ("  -v  : Write to Logfile [off]\n");
  print ("  -B  : Save blast output (in $tmpdir) [off]\n");
  print ("  -b  : Show alignments for [$alignments]\n");
  print ("  -R  : Store raw blast profile [off]\n");
  print ("  -a  : Blast database chunks [$chunks]\n");
  print ("  -o  : output profile name\n");
  print ("\n");
 exit;
} # Usage

#
# Open input
# If no file name is given, use standard input
#
if (not defined($Getopt::Std::opt_i)){
  # Read from standard input
  *INP = *STDIN;
} 
else{
  # Read from file
  if (($Getopt::Std::opt_i=~/\.gz$/) || ($Getopt::Std::opt_i=~/\.Z$/)){
    open(INP,"gunzip -c $Getopt::Std::opt_i |") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
  else{
    open(INP,"<$Getopt::Std::opt_i") || die ("can't open file $Getopt::Std::opt_i: $!");
  }
}

if (defined($Getopt::Std::opt_v)){
  if (defined($Getopt::Std::opt_L)){
    open(LOGFILE,">$Getopt::Std::opt_L");
  }
  else{
    *LOGFILE = *STDERR;
  }
}

#
# Blast database
#
if (defined($Getopt::Std::opt_d)){
  $db = $Getopt::Std::opt_d;
}

if (defined($Getopt::Std::opt_j)){
  $itr=$Getopt::Std::opt_j;
}

if (defined($Getopt::Std::opt_e)){
  $evalue=$Getopt::Std::opt_e;
}
if (defined($Getopt::Std::opt_t)){
  $tmpdir = $Getopt::Std::opt_t;
}
if (! -d "$tmpdir"){
  my $cmd = "mkdir -m 0777 -p $tmpdir";
  if (defined($Getopt::Std::opt_v)){
    print LOGFILE ("# $cmd\n");
  }
  system("$cmd");
}
if (defined($Getopt::Std::opt_b)){
    $alignments=$Getopt::Std::opt_b;
}
if (defined($Getopt::Std::opt_a)){
  $chunks=$Getopt::Std::opt_a;
}
###################################################################
my $orig_dir = cwd();
if (defined($Getopt::Std::opt_v)){
  print LOGFILE ("# Working Directory:\t$orig_dir\n");
  print LOGFILE ("# Command:\t\t$command\n");
}

if (defined($Getopt::Std::opt_I)){
    $TYPE=$Getopt::Std::opt_I;
}

#####################################################################
#
# Main programme
#
#####################################################################
my $readTYPE = "read".$TYPE;
my $i=0;
my @names=();
my @allnames=();
my %r={};
my %accnames={};
while (! eof (INP)){
  my %rec = {};

  #
  # Read input file
  #
  %rec = &$readTYPE(\*INP);
  
  if ($TYPE eq "FASTA"){
      for (my $i=1;$i<=$rec{len};$i++){
	  $rec{sec}[$i]='.';
      }
  }
  my $Qname = $rec{name};
  if (exists($accnames{$Qname})){
      print ("\nA duplicate sequence identifier was used: '$Qname'.\nOnly unique sequence identifiers are allowed!\n");
	exit 1;
  }
  $accnames{$Qname}=1;
  push(@allnames, $Qname);
  my $file="$Qname".".prof";

  #
  # Make the fasta file
  #
  $i++;
  my $fasta = "$tmpdir/$Qname.fasta";
  my $how = "$tmpdir/$Qname.how";

  push(@names,$Qname);
  open(FASTA, ">$fasta");
  writeFASTA(\%rec, \*FASTA);
  close(FASTA);

  open(HOW, ">$how");
  writeHOW(\%rec, \*HOW);
  close(HOW);
  print LOGFILE ("# $i :\tMade $fasta and $how\n");
}
close(INP);
my $entries=$i;
my $acids=20;

my $k=0;
my $q_name="";

foreach $q_name (@names){
  $k++;
  #
  # Run blastpgp
  #
  my $raw = "$profile_db_lookup/$q_name.raw";
  my $fasta = "$tmpdir/$q_name.fasta";
  my $blastout = "$tmpdir/$q_name.blastout.gz";
  my $blast_profile = "$tmpdir/$q_name.raw";
  if (-e "$raw"){
    system("cp $raw $blast_profile");
  }

  $cmd = "$blastpgp  -d $db -e $evalue -a $chunks -j $itr -b $alignments -Q $blast_profile -i $fasta> /dev/null";

  if (defined($Getopt::Std::opt_B)){
      $cmd = "$blastpgp -d $db -e $evalue -a $chunks-j $itr -b $alignments -Q $blast_profile -i $fasta | gzip > $blastout";
  }

  system("$cmd");
}

my $i=0;

foreach $q_name (@names){
  $i++;
  my %rec={};
  my %prof = {};

  my $blast_profile = "$tmpdir/$q_name.raw";
  my $how = "$tmpdir/$q_name.how";

  open(HOW,"<$how");
  %rec=readHOW(\*HOW);
  close(HOW);

  if (-e "$blast_profile"){
      open(PROF, "<$blast_profile");
      %prof = readBLASTPROF(\*PROF);
#    $prof{name} = substr($q_name, 0, 13);
      $prof{name} = $q_name;
      for (my $j = 1; $j <= $prof{len}; $j++){
	  if (defined($rec{sec}[$j])){
	      $prof{sec}[$j] = $rec{sec}[$j];
	  }
	  else{
	      $prof{sec}[$j] = ".";
	  }
	  if (($prof{sec}[$j] eq '') || ($prof{sec}[$j] eq ' ')){
	      $prof{sec}[$j] = ".";
	  }
      }
    open(PROFDB,">$tmpdir/$q_name.prof");
    writePROF(\%prof, $acids, \*PROFDB);
    close(PROFDB);
    $r{$q_name}="cat $tmpdir/$q_name.prof";
  }
  else{
    #
    # Profile will be made via a blosum62 matrix
    #
    %prof = makeBLOSUM62PROF(\%rec);
    open(PROFDB,">$tmpdir/$q_name.prof");
    writePROF(\%prof, $acids, \*PROFDB);
    close(PROFDB);
    $r{$q_name}="cat $tmpdir/$q_name.prof";
  }
}
$i=0;
foreach my $id (@allnames){
  $i++;
  my $direct = ">>";
  if ($i == 1){
    $direct = ">";
  }
  my $cmd = "$r{$id}";
  if (defined($Getopt::Std::opt_o)){
    $cmd .= " $direct $Getopt::Std::opt_o";
  }
  if (defined($Getopt::Std::opt_v)){
    print LOGFILE ("$cmd\n");
  }
  system("$cmd");
}
#
# Cleanup
#
if ( (-d "$tmpdir") && (! defined($Getopt::Std::opt_k))){
  my $cmd = "rm -r $tmpdir";
  if (defined($Getopt::Std::opt_v)){
    print LOGFILE ("# Doing '$cmd'\n");
  }
  system("$cmd");
}


