package ACT;
use Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(readACT writeACT write_noneACT statACT mcc);
use strict;
#######################################################################
=pod

=head2 readACT

        readACT - read in ACT format

=head3 DESCRIPTION

        readACT - read a profile entry

=head3 SYNOPSIS

        %rec=readACT(\*filehandle)

=head3 PARAMETERS

        \*filehandle              Handle to file or STDOUT
	$rec{len}                 length of profile
	$rec{name}                name of profile
	$rec{sec}[$i]             secondary structure
	$rec{seq}[$i]             primary sequence
	$rec{resnum}[$i]          residue number $i
	$rec{act}[$i][$k]         Activity for element $k
=cut#######################################################################
sub readACT{
    my ($fh) = @_;
    my %rec={};
    my $i=0;
    my $first="true";
    my $name="not_important_now";
    my $previous_name="not_important_now";
    my $pos=0;
    while (($name eq $previous_name) && (!eof($fh)))
      {
	$_=<$fh>;
	if (/^#/){
	    next;
        }
	$i++;
	chomp($_);
	my @words=();
	@words=split(/\s+/,$_);
	$name=$words[2];
	
	if ($first eq "true"){
	  $previous_name=$name;
	  $first="false";
	}
	
	if ("$name" eq "$previous_name"){
	  $previous_name=$name;
	  $pos=tell $fh;
	}
	else{
	  seek $fh,$pos,0;
	  return %rec;
	}

	$rec{sec}[$i]=$words[0];
	$rec{seq}[$i]=$words[1];
	$rec{name}=$words[2];
	$rec{resnum}[$i]=$words[3];

	my $k=1;
	while (defined($words[3+$k])){
	  $rec{act}[$i][$k]=$words[3+$k];
	  $k++;
	}
	$rec{len}=$i;
	$previous_name=$rec{name};
      }
  return %rec;
}
######################################################################
=pod

=head2 writeACT

        writeACT - writes a profile

=head3 DESCRIPTION

        writeACT - writes a profiles

=head3 SYNOPSIS

        writeACT (\%rec,$columns, \*filehandle)

=head3 PARAMETERS

        \*filehandle    Handle to file or STDOUT
        \%rec         Reference to a hash table.
        $columns        Is the number of columns to write   
=cut#######################################################################
sub writeACT{
  my ($rh_rec, $fh,$format2) = @_;

  for (my $i=1;$i<=$rh_rec->{len};$i++){
    if (! defined($rh_rec->{resnum}[$i])){
	$rh_rec->{resnum}[$i]=$i;
      }
      printf $fh ("%1s %1s  %-18s %5d ",$rh_rec->{sec}[$i],$rh_rec->{seq}[$i],$rh_rec->{name},$rh_rec->{resnum}[$i]);
      my $k=1;
        while (defined($rh_rec->{act}[$i][$k])){
	  printf $fh (" %7.3f",$rh_rec->{act}[$i][$k]);
	  $k++;
        }
        printf $fh ("\n");
    }
  return;
  }


sub write_noneACT{
  my ($rh_rec, $fh) = @_;
  
  for (my $i=1;$i<=$rh_rec->{len};$i++){
    if (! defined($rh_rec->{resnum}[$i])){
      $rh_rec->{resnum}[$i]=$i;
    }
    printf $fh ("%1s %1s  %-18s %5d\n",$rh_rec->{sec}[$i],$rh_rec->{seq}[$i],$rh_rec->{name},$rh_rec->{resnum}[$i]);
  }
  return;
}


sub statACT{
  my ($rh_rec, $arr) = @_;
  my %stat={};
  my $cat="";
  my $skip_cat = "-";

  foreach $cat (@{$arr}){
    $stat{tp}{$cat}=0;
    $stat{tn}{$cat}=0;
    $stat{fp}{$cat}=0;
    $stat{fn}{$cat}=0;
    $stat{total}{$cat}=0;
  }
  my $cat_num=-1;

  for (my $i=1;$i<=$rh_rec->{len};$i++){
    my $k=1;
    my $max=-999;
    if ($rh_rec->{sec}[$i] eq "-"){
      next;
    }
    while (defined($rh_rec->{act}[$i][$k])){
      if ($rh_rec->{act}[$i][$k] > $max){
	$max=$rh_rec->{act}[$i][$k];
	$cat_num=$k -1;
      }
      $k++;
    }
    $stat{total}{$rh_rec->{sec}[$i]}++;
    #print ("$rh_rec->{sec}[$i]\t$stat{total}{$rh_rec->{sec}[$i]}\n");
    if ($$arr[$cat_num] eq $rh_rec->{sec}[$i]){
      foreach $cat (@{$arr}){
	if ($cat ne $$arr[$cat_num]){
	  $stat{tn}{$cat}++;
	}
	else{
	  $stat{tp}{$cat}++;
	}
      }
    }
    else{
      foreach $cat (@{$arr}){
	if ($cat ne $$arr[$cat_num]){
	  $stat{fn}{$cat}++;
	}
	else{
	  $stat{fp}{$cat}++;
	}
      }
    }
  }
  return(%stat);
}

sub mcc{
  my ($tp, $tn, $fp, $fn) = @_;
  my $corr=0;
  my $perc=0;
  my $a = $tn + $fn;
  my $b = $tn + $fp;
  my $c = $tp + $fn;
  my $d = $tp + $fp;

  my $lower =  sqrt($a * $b * $c * $d);
  $corr = ($tp * $tn - $fn* $fp)/$lower;
  my $test = $tp + $fn;
  if ($test == 0){
    $perc=100.0;
  }
  else{
    $perc=$tp *100 /($tp + $fn);
  }
  return($corr,$perc);
}
1;
