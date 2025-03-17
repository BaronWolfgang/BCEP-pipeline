package PROF;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(readPROF writePROF makeBLOSUM62PROF blosum62);
use strict;
#######################################################################
=pod

=head2 readPROF

        readPROF - read in PROF format

=head3 DESCRIPTION

        readPROF - read a profile entry

=head3 SYNOPSIS

        %rec=readPROF(\*filehandle)

=head3 PARAMETERS

        \*filehandle              Handle to file or STDOUT
	$rec{len}                 length of profile
	$rec{name}                name of profile
	$rec{sec}[$i]             secondary structure
	$rec{seq}[$i]             primary sequence
	$rec{resnum}[$i]          residue number
	$rec{score}[$i][$k]       profile score for element $k 
	$rec{activity}[$i][$k]    Activity for element $k
=cut#######################################################################
sub readPROF{
    my ($fh) = @_;
    my %rec;
    my $i=0;
    my $first="true";
    my $counter=-999;
    my $counter_previous=-9999;
    my $name="not_important_now";
    my $previous_name="not_important_now";
    my $pos=0;
    while (!eof($fh))
      {
	$i++;
	my @words=split(/\s+/,<$fh>);
	$name=$words[2];
	$counter=$words[3];
	
	if ($first eq "true")
	  {
	    $previous_name=$name;
	    $first="false";
	    $counter_previous=-9999;
	  }
	if (($name eq $previous_name) && ($counter > $counter_previous))
	  {
	    #$previous_name=$name;
	    $counter_previous=$counter;
	    $pos=tell $fh;
	  }
	else
	  {
	    seek $fh,$pos,0;
	    return %rec;
	  }
	$rec{sec}[$i]=$words[0];
	$rec{seq}[$i]=$words[1];
	$rec{name}=$words[2];
	$rec{resnum}[$i]=$words[3];
	for (my $k=1;$k<=20;$k++)
	  {
	    $rec{score}[$i][$k]=$words[3+$k];
	  }
	my $k=1;
	while (defined($words[23+$k]))
	{
	  $rec{act}[$i][$k]=$words[23+$k];
	  $k++;
	}
	$rec{len}=$i;
#	$previous_name=$rec{name};
	$previous_name=$name;
#	$counter_previous=$counter;
      }
  return %rec;
}
######################################################################
=pod

=head2 writePROF

        writePROF - writes a profile

=head3 DESCRIPTION

        writePROF - writes a profiles

=head3 SYNOPSIS

        writePROF (\%rec,$columns, \*filehandle)

=head3 PARAMETERS

        \*filehandle    Handle to file or STDOUT
        \%rec         Reference to a hash table.
        $columns        Is the number of columns to write   
=cut#######################################################################
sub writePROF{
  my ($rh_rec,$columns,$fh) = @_;

    for (my $i=1;$i<=$rh_rec->{len};$i++){
      if (! defined($rh_rec->{resnum}[$i])){
	$rh_rec->{resnum}[$i]=$i;
      }
      printf $fh ("%1s %1s %-20s %6d ",$rh_rec->{sec}[$i],$rh_rec->{seq}[$i],$rh_rec->{name},$rh_rec->{resnum}[$i]);
      for (my $k=1;$k<=$columns;$k++){
	printf $fh ("%3d ",$rh_rec->{score}[$i][$k]);
      }

#      $k=1;
#      while (defined($rh_rec->{score}[$i][$k])){
#	printf $fh ("%3d ",$rh_rec->{score}[$i][$k]);
#	$k++;
#      }

      my $k=1;
      while (defined($rh_rec->{act}[$i][$k])){
	printf $fh (" %6.4f",$rh_rec->{act}[$i][$k]);
	$k++;
      }
      printf $fh ("\n");
    }
  return;
  }
#############################################################################################
sub num2aa{
  my $i = 0;
  my @vec=();
  my @aa = qw(A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V X);

  foreach my $id (@aa){
    $i++;
    $vec[$i] = $id;
  }
  return(@vec);
}


sub aa2num{
  my $i = 0;
  my %aa;
  $aa{A} = 1;
  $aa{R} = 2;
  $aa{N} = 3;
  $aa{D} = 4;
  $aa{C} = 5;
  $aa{Q} = 6;
  $aa{E} = 7;
  $aa{G} = 8;
  $aa{H} = 9;
  $aa{I} = 10;
  $aa{L} = 11;
  $aa{K} = 12;
  $aa{M} = 13;
  $aa{F} = 14;
  $aa{P} = 15;
  $aa{S} = 16;
  $aa{T} = 17;
  $aa{W} = 18;
  $aa{Y} = 19;
  $aa{V} = 20;

  return(%aa);
}

sub makeBLOSUM62PROF {
  my ($rh_rec) = @_;
   
   my %blosum = ();
   my %prof;
 
   $prof{name} = $rh_rec->{name};
   $prof{len} = $rh_rec->{len};
   for (my $i = 1; $i <= $rh_rec->{len}; $i++){
     my $aa = $rh_rec->{seq}[$i];
     %blosum = blosum62($aa);
   
     $prof{seq}[$i] = $rh_rec->{seq}[$i];
     $prof{resnum}[$i] = $rh_rec->{resnum}[$i];

     if (defined($rh_rec->{sec}[$i])){
       $prof{sec}[$i] = $rh_rec->{sec}[$i];
     }
     else{
       $prof{sec}[$i] = ".";
     }

     for (my $j = 1; $j <= 20; $j++){
       $prof{score}[$i][$j] = $blosum{$aa}[$j];
     }
     #print ("\n");
   }
  return(%prof);
}

sub blosum62{
  my ($acid) = @_;
  my %vec;
  #
  #            A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
  #
  @{$vec{A}} = qw(4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4);
  @{$vec{R}} = qw(-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4);
  @{$vec{N}} = qw(-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4);
  @{$vec{D}} = qw(-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4);
  @{$vec{C}} = qw( 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4);
  @{$vec{Q}} = qw(-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4);
  @{$vec{E}} = qw(-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4);
  @{$vec{G}} = qw(0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4);
  @{$vec{H}} = qw(-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4);
  @{$vec{I}} = qw(-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4);
  @{$vec{L}} = qw(-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4);
  @{$vec{K}} = qw(-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4);
  @{$vec{M}} = qw(-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4);
  @{$vec{F}} = qw(-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4);
  @{$vec{P}} = qw(-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4);
  @{$vec{S}} = qw(1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4);
  @{$vec{T}} = qw(0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4);
  @{$vec{W}} = qw(-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4);
  @{$vec{Y}} = qw(-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4);
  @{$vec{V}} = qw(0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4);
  @{$vec{B}} = qw(-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4);
  @{$vec{Z}} = qw(-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4);
  @{$vec{X}} = qw(0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4);
  @{$vec{"*"}} = qw(-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 1);

  if ($acid=~/[ARNDCQEGHILKMFPSTWYVBZX\*]/){
    #print ("$acid\t'$vec{$acid}[0]'\n");
    return(%vec);
#    return($vec{$acid});
  }
  else{
    print STDERR ("# Error in blosum62 PROF.pm: unknown amino acid '$acid'\n");
    return;
  }
}

1;
