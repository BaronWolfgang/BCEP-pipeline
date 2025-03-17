#!/usr/bin/perl
use strict;
use Getopt::Std;
use lib "$ENV{NetSurfP}/modules";
use PROF;
use ACT;
my $winsize=13;
my $neurons=20;
#
# Process command line
#
getopts('hi:o:w:t:n:j:v')||Usage();
#
# Usage
#
if (defined($Getopt::Std::opt_h)||defined($Getopt::Std::opt_h)){
  # Print help message
  Usage();
}

sub Usage {
  print ("Usage: $0 [-h] [-i file] [-o file] [-w number]\n");
  print ("Description:\n");
  print ("number2win.pl - generate windows for howlin given a profile\n");
  print ("\n");
  print ("Options:\n");
  print ("  -h  : display this message\n");
  print ("  -i  : input profile file\n");
  print ("  -j  : input act file\n");
  print ("  -n  : neurons per character or line [$neurons]\n");
  print ("  -w  : window size [$winsize]\n");
  print ("  -o  : output howlin file [STDOUT]\n");
  print ("  -v  : Verbose mode [off]\n");
  print ("\n");
 exit;
}

#
# Open input profile
#
if (not defined($Getopt::Std::opt_i)){
    print STDERR ("# Option -i has not been defined\n");
    exit
} 
else{
    # Read from file
    if (($Getopt::Std::opt_i=~/\.gz$/) || ($Getopt::Std::opt_i=~/\.Z$/)){
	open(PROF,"gunzip -c $Getopt::Std::opt_i |") || die ("can't open file $Getopt::Std::opt_i: $!");
    }
    else{
	open(PROF,"<$Getopt::Std::opt_i") || die ("can't open file $Getopt::Std::opt_i: $!");
    }
}
#
# Open input activity file
#
if (not defined($Getopt::Std::opt_j)){
    print STDERR ("# Option -j has not been defined\n");
    exit
} 
else{
    # Read from file
    if (($Getopt::Std::opt_j=~/\.gz$/) || ($Getopt::Std::opt_j=~/\.Z$/)){
	open(ACT,"gunzip -c $Getopt::Std::opt_j |") || die ("can't open file $Getopt::Std::opt_j: $!");
    }
    else{
	open(ACT,"<$Getopt::Std::opt_j") || die ("can't open file $Getopt::Std::opt_j: $!");
    }
}
#
# If not file name is given, use standard output
#
if (not defined($Getopt::Std::opt_o)){
  # Output goes to std output
  *OUT = *STDOUT;
} else {
  # Open file to write to
  open(OUT, ">$Getopt::Std::opt_o");
}
if (defined($Getopt::Std::opt_w)){
    $winsize=$Getopt::Std::opt_w;
}
###############################################################################
# Main
#
###############################################################################
my $target_value=1;
my $val;
my $last_neuron=9999;
while (! eof (PROF)){
  my %rec=();
  my %activ=();
  my $col=20;
  %rec=readPROF(\*PROF);
  %activ=readACT(\*ACT);
  if (defined($Getopt::Std::opt_v)){
    print STDERR ("# Reading new seq from profile: $rec{name}\n");
  }
  my $half=int($winsize/2);
  for (my $i=1;$i<=$rec{len};$i++){
      for (my $k=$i-$half;$k<=$i+$half;$k++){
	  for (my $j=1;$j<=$neurons;$j++){
	      if (($k<=0) || ($k>$rec{len})){		  
		  $val=0.0;
	      }
	      else{
		  #
		  # Profile scores are multiplied with a factor.
		  # The factor has previously been found to be the optimum
		  # value for secondary structure prediction
		  #
		  $val=$rec{score}[$k][$j]*0.075;
	      }
	      printf OUT ("%6.3f ",$val);
	  }
	  for (my $j=1;$j<=2;$j++){
	      if (($k<=0) || ($k>$rec{len})){		  
		  $last_neuron=1.0;
		  $val=0.0;
	      }
	      else{
		  $val = $activ{act}[$k][$j];
		  $last_neuron=0.0;
	      }
	      printf OUT ("%6.3f ",$val);
#	      print STDERR ("\n# k=$k\tj=$j\tval='$val'\n");
	  } 
	  printf OUT ("%6.3f ",$last_neuron);
     }    
      printf OUT ("%6.4f\n",$target_value);
  }
}


#############################################################################
