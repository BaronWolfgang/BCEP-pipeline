#######################################################################
=pod

=head1 NAME

        FASTA - Support of FASTA format.

=cut
#######################################################################

package FASTA;                     # Package name              

use Exporter;                       # Load Exporter module (APP, p91)
@ISA = qw(Exporter);                # Inherit from Exporter module
@EXPORT = qw(readFASTA writeFASTA remove_X_FASTA npFASTA writeBLATFORMAT); # List of symbols, OK to export
use strict;
#######################################################################
=pod

=head2 writeFASTA

        writeFASTA - write in FASTA format

=head3 DESCRIPTION

        writeFASTA writes a protein sequence in the FASTA format.

=head3 SYNOPSIS

        writeFASTA(\%FASTAentry, \*filehandle)

=head3 PARAMETERS

        \%rec		  Reference to FASTA hash
        \*fh		  A reference to filehandle (file or STDOUT)

         eg writeFASTA(\%rec, \*fh)

        The following fields must be set in the calling program

        $rec{len}   length
        $rec{name}   name of sequence
        $rec{desc}   description of sequence
        $rec{seq}[$i]    The sequence

=cut
#######################################################################

sub writeFASTA {
  my ($rh_rec, $fh)  = @_;

  #
  # WRITE THE FASTA HEADER
  #

  printf $fh (">%s", $rh_rec->{name});
  if (defined($rh_rec->{desc}) && ($rh_rec->{desc} ne "(null)")){
    printf $fh (" %s", $rh_rec->{desc});
  }
  printf $fh ("\n");

  #
  # WRITE THE FASTA SEQUENCE
  #
  my $i = 1;
  
  while (defined($rh_rec->{len}) && ($i<=$rh_rec->{len})){
    my $acid="";

    if (defined($rh_rec->{seq}[$i])){
      $acid =  $rh_rec->{seq}[$i];
#      print STDERR ("# $i\t$acid\n");
    }	
    else{
      print STDERR ("# Undefined residue '$acid' no $i in $rh_rec->{name} (FASTA.pm:writeFASTA)\n");
      $rh_rec->{len}--;
      next;
    }
    printf $fh ("%s", $acid);
    if (($i % 60) == 0) {
      print $fh ("\n");
    }
    $i++;
  }
    
  if ((int($rh_rec->{len}) % 60) != 0) {
    print $fh ("\n");
  }
 undef($rh_rec);
} # writeFASTA

##################################################
=pod

=head2 readFASTA

        readFASTA - read a FASTA entry from an open file

=head3 DESCRIPTION

        readFASTA reads a FASTA entry

=head3 SYNOPSIS

        %rec = readFASTA(\*fh)

=head3 PARAMETERS

        %rec     	 A FASTA hash
        \*fh     	 A reference to a filehandle (file or STDIN)

        $rec{len}    	length
        $rec{name}   	name of sequence
        $rec{desc}   	description of sequence
        $rec{seq}[$i]   sequence

=cut
##################################################

sub readFASTA {
  my ($fh) = @_;
  my %rec = ();
  my @newsequence = ();
  my @sequence = ();

  $_ = <$fh>;
  #
  # Skip comment lines
  #
  while (/^\#/){
    $_ = <$fh>;
    #chomp($_);
  }

  #
  # Skip lines not starting with ">"
  #
  while (substr($_,0,1) ne ">"){
    $_ = <$fh>;
  }

  if (substr($_, 0, 1) ne ">" )
    {
      print STDERR ("Incorrect FASTA format\n");
      $rec{error}=1;
      return %rec;
      exit;
    }

  while (substr($_,0,1) eq ">"||substr($_,0,1) eq " ") {
    $_ = substr($_,1);
    $_=~ s/\|/_/g;
  }
  
  #
  # THE FIRST WORD IN THE COMMENT LINE BECOMES THE NAME OF THE ENTRY
  #
  $_=~ s/>/_/g;

  my @fields=split(/\s+/,$_);
  $rec{name}=$fields[0];
  if ($rec{name} eq ''){

	  $rec{name} = "Sequence";
  };
  
  #
  # THE REMAINING PART OF THE COMMENT LINE BECOMES ADDED
  # 

  my $start=length($rec{name});
  my $end=length($_);
  $rec{desc}=substr($_, $start, $end);
  chomp($rec{desc});


  #
  # Read the sequence
  #
  $_ = <$fh>;
  chomp($_);
  $_=~ s/\W//g;
  $_=~ s/\s+//g;
  $_=~ s/\r//g;
  $_=~ s/\n//g;

  #
  # Change the sequence to uppercase letters
  #
  my $tmpstring = uc($_);
  $_=$tmpstring;
  my $len=0;
  #
  # Check if the previous fasta entry was empty
  #
  if (substr($_, 0, 1) eq ">")
  {
    print STDERR "#\nTHE SEQUENCE: $rec{name} CONTAINS NO RESIDUES\n#\n";
    $rec{error}=1;
    return %rec;
    exit;
  }

  @sequence = split(//, "");
  $rec{seq}=\@sequence;
  $sequence[0]=" ";
  my $newentry=0;
  my $endoffile=0;
  $rec{len}=0;
  while ((! $newentry) && (! $endoffile)) 
  {
    chomp($_);
    @newsequence = split(//,$_);
    
    push (@sequence, @newsequence);
    $rec{len}=$rec{len}+length($_);
    
    if (! eof($fh))
    {
      $_ = <$fh>;
      $len=length($_);
      if (substr($_,0,1) eq ">")
      {
        $len=length($_);
        $newentry=1;
      }
    }
    else
    {
      $endoffile=1;
    }
  }

  if ($newentry)
  {
    $len=$len*-1;
    seek $fh,$len,1;
  }
  return %rec;
} # readFASTA

sub remove_X_FASTA{
  my ($rh_rec) = @_;

  my %new = {};

  $new{name} = $rh_rec->{name};
  $new{desc} = $rh_rec->{desc};
  my $k=0;
  for (my $i=1; $i<=$rh_rec->{len}; $i++){
    if( ($rh_rec->{seq}[$i] ne "X") && ($rh_rec->{seq}[$i] ne "x") ){
      $k++;
      $new{seq}[$k] = $rh_rec->{seq}[$i];
      $new{len}++;
    }
  }
  return (%new),
}

###############################################################################################
#
# NUCLEOTIDE CONVERSION TABLE
#
###############################################################################################
sub npFASTA{
  my (%np) = @_;
$np{A}{A}{A} = 'K';
$np{A}{A}{T} = 'N';
$np{A}{A}{C} = 'N';
$np{A}{A}{G} = 'K';
$np{A}{T}{A} = 'I';
$np{A}{T}{T} = 'I';
$np{A}{T}{C} = 'I';
$np{A}{T}{G} = 'M';
$np{A}{C}{A} = 'T';
$np{A}{C}{T} = 'T';
$np{A}{C}{C} = 'T';
$np{A}{C}{G} = 'T';
$np{A}{G}{A} = 'R';
$np{A}{G}{T} = 'S';
$np{A}{G}{C} = 'S';
$np{A}{G}{G} = 'R';


#
# First is T
#
$np{T}{A}{A} = '*';
$np{T}{A}{T} = 'Y';
$np{T}{A}{C} = 'Y';
$np{T}{A}{G} = '*';
$np{T}{T}{A} = 'L';
$np{T}{T}{T} = 'F';
$np{T}{T}{C} = 'F';
$np{T}{T}{G} = 'L';
$np{T}{C}{A} = 'S';
$np{T}{C}{T} = 'S';
$np{T}{C}{C} = 'S';
$np{T}{C}{G} = 'S';
$np{T}{G}{A} = '*';
$np{T}{G}{T} = 'C';
$np{T}{G}{C} = 'C';
$np{T}{G}{G} = 'W';


#
# First is C
#
$np{C}{A}{A} = 'Q';
$np{C}{A}{T} = 'H';
$np{C}{A}{C} = 'H';
$np{C}{A}{G} = 'Q';
$np{C}{T}{A} = 'L';
$np{C}{T}{T} = 'L';
$np{C}{T}{C} = 'L';
$np{C}{T}{G} = 'L';
$np{C}{C}{A} = 'P';
$np{C}{C}{T} = 'P';
$np{C}{C}{C} = 'P';
$np{C}{C}{G} = 'P';
$np{C}{G}{A} = 'R';
$np{C}{G}{T} = 'R';
$np{C}{G}{C} = 'R';
$np{C}{G}{G} = 'R';



#
# First is C
#
$np{G}{A}{A} = 'E';
$np{G}{A}{T} = 'D';
$np{G}{A}{C} = 'D';
$np{G}{A}{G} = 'E';
$np{G}{T}{A} = 'V';
$np{G}{T}{T} = 'V';
$np{G}{T}{C} = 'V';
$np{G}{T}{G} = 'V';
$np{G}{C}{A} = 'A';
$np{G}{C}{T} = 'A';
$np{G}{C}{C} = 'A';
$np{G}{C}{G} = 'A';
$np{G}{G}{A} = 'G';
$np{G}{G}{T} = 'G';
$np{G}{G}{C} = 'G';
$np{G}{G}{G} = 'G';
  return(%np);
}

#####################################################################

sub writeBLATFORMAT {
  my ($rh_rec, $start, $mp, $fh)  = @_;

  #
  # $mp is minus or plus direction
  #
  
  #
  # WRITE THE SEQUENCE
  #
  my $i=0;
  if ($mp eq "+"){
    $i = 1;
  }
  else{
    $i=-1;
  }
  my $k = 0;
  my $ii=1;
  my $number=0;
  my $acid="";
  while (defined({$rh_rec->{len}}) & ($ii<={$rh_rec->{len}}))
    {
      $number = $start + $i;
      if (defined({$rh_rec->{seq}[$ii]})){
	$acid =  %{$rh_rec->{seq}[$ii]};
      }	
      else{
       print STDERR ("# Undefined residue '$acid' no $ii in %{$rh_rec->{name}}\n");
       exit;
      }
      printf $fh ("%s", $acid);

      if (($ii % 10) == 0) {
        print $fh (" ");
      }

      if (($ii % 50) == 0) {
	$k = 0;
        print $fh (" $number\n");
      }
      if ($mp eq "+"){
	$i++;
      }
      else{
	$i--;
      }
      $k++;
      $ii++;
    }


  my $j=0;
  if ($k != 0){
    $k--;
  }

  if ((int(%{$rh_rec->{len}}) % 50) != 0) {
    for ($i=$number+$j; $k<=50; $k++){
      $j++;
      if (($i % 10) == 0) {
        print $fh ("  ");
      }
      else{
	print $fh (" ");
      }
    }
    print $fh (" $number\n");
  }
 undef ($rh_rec);
} # writeBLATFORMAT
1;

##################################################
