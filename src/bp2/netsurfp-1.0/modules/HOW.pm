package HOW;

use Exporter;
@ISA = qw(Exporter);        
@EXPORT = qw(writeHOW readHOW);
use strict;

sub writeHOW {
  my ($rh_rec, $fh)  = @_;
  # write the header
  printf $fh ("%6d %s", $rh_rec->{len}, $rh_rec->{name});
  if (defined($rh_rec->{desc})){
    chomp($rh_rec->{desc});
    printf $fh (" %s", $rh_rec->{desc});
  }
  printf $fh ("\n");

  # write the sequence
  my $i = 1;
  while (defined($rh_rec->{len}) && ($i<=$rh_rec->{len}))
    {
      if (defined($rh_rec->{seq}[$i])){
        printf $fh ("%s", $rh_rec->{seq}[$i]);
       }
      if (($i % 80) == 0){
        print $fh ("\n");
      }
      $i++;
    } 
  if ((int($rh_rec->{len}) % 80) != 0){
    print $fh ("\n");
  }

  # write the secondary structure
  $i=1;
  while ($i <= $rh_rec->{len}) {
    if ($rh_rec->{seq}[$i] eq "-") {
      printf $fh ("-");
    } elsif ($rh_rec->{sec}[$i] eq "none") {
      printf $fh (".");
    } else {
      printf $fh ("%s", $rh_rec->{sec}[$i]);
    }
    if (($i % 80) == 0){
      print $fh ("\n");
    }
    $i++;
  }
  if ((int($rh_rec->{len}) % 80) != 0) {
    print $fh ("\n");
  }
} 


sub readHOW {
  my ($fh) = @_;

  my %rec=();
  my @fields=();
  my $line;
  my $i;
  my $j;
  my $nolines;
  my $rest;

  $line = <$fh>;
  chomp($line);
  @fields=split/\s+/,$line;
  if ($fields[0] eq ""){
    $rec{len}=$fields[1];
    $rec{name}=$fields[2];
#    $rec{name}=~ tr/a-z/A-Z/;
#    print STDERR ("# '$rec{name}'\n");
#    print STDERR ("# '$rec{len}'\n");
#    print STDERR ("# '$line'\n");
    $rec{desc}="";
    for ($i=3;$i<scalar(@fields);$i++){
      $rec{desc} .= " ".$fields[$i];
      chomp($fields[$i]);
#      print STDERR ("'$fields[$i]\n'");
    }
  }
  else
  {
    $rec{len}=$fields[0];
    $rec{name}=$fields[1];
    $rec{desc}="";
    for ($i=2;$i<scalar(@fields);$i++){
      chomp($fields[$i]);
      $rec{desc} .= "$fields[$i] ";
    }
    chop($rec{desc});
  }
  $nolines=int((int($rec{len})/80));
  # read sequence
  $i=1;
  while ($i<=$nolines) {
    $line=<$fh>;
    chomp($line);
    $j=0;
    while ($j<80) {
      $rec{seq}[$j+1+80*($i-1)]=substr($line, $j, 1);
      $j++;
    }
    $i++;
  }
  # read last line (shorter)
  $rest=int($rec{len})-80*int($nolines);
  if ($rest > 0) {
    $line=<$fh>;
    chomp($line);
  }
  $i=0;
  while ($i<$rest) {
    $rec{seq}[int($i)+1+80*int($nolines)]=substr($line, $i, 1);
    $i++;
  }

  # read secondary structure
  $i=1;
  while ($i<=$nolines) {
    $line=<$fh>;
    chomp($line);
    $j=0;
    while ($j<80) {  # problem: last line is shorter
      $rec{sec}[$j+1+80*($i-1)]=substr($line, $j, 1);
      $j++;
    }
    $i++;
  }
  # read last line (shorter)
  $rest=int($rec{len})-80*int($nolines);
  if ($rest > 0) {
    $line=<$fh>;
    chomp($line);
  }
  $i=0;
  while ($i<$rest) {
    $rec{sec}[$i+1+80*int($nolines)]=substr($line, $i, 1);
    $i++;
  }

  return %rec;
} 
1;
