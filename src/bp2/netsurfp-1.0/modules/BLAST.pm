#######################################################################
=pod

=head1 NAME

        BLAST - Support of BLAST format.

=cut
#######################################################################

package BLAST;                     # Package name              

use Exporter;                       # Load Exporter module (APP, p91)
@ISA = qw(Exporter);                # Inherit from Exporter module
@EXPORT = qw(readBLAST writeBLAST_short writeBLAST_long writeBLAST readBLASTPROF writeBLAST_comment); # List of symbols, OK to export

#######################################################################
=pod

=head2 readBLAST

        readBLAST - read a BLAST entry from an open file

=head3 DESCRIPTION

        readBLAST reads a BLAST entry

=head3 SYNOPSIS

        %rec = readBLAST(\*fh)

=head3 PARAMETERS

        %rec 		BLAST entry
        \*fh      	File handle  to file or STDIN

        $rec{len}    	length of query sequence
        $rec{name}   	name of query sequence
	a.o.
=cut
##################################################

sub readBLAST {
#  my ($fh, $idx_file, $replace) = @_;
  my ($fh, $main_hits) = @_;
  my %rec;
  
  #
  # Initialize entries
  #
  $rec{name} = "-";
  $rec{Qhits} = 0;
  $rec{Qlen} = 0;
  $rec{Qhits} = 0;
  $i = 0;
  $MAXLEN = 200;
  $found = 0;
  $type = "Protein";

  #
  # Parse the file
  #
  $_= <$fh>;
  while ( (! eof ($fh)) && !(/^S2:/) && ($rec{Qhits} < $main_hits)){
    print ("$_'$rec{Qhits}'\t'$main_hits'\n");
    if (eof($fh)){
      print ("## returning\n");
      return(%rec);
    }
    $_= <$fh>;
    if ($found == 1){
       $line++;
       $rec{line}[$line] = $_;
       $rec{lines}++;
    }
    
    if (/Query= (\S+)/){
      $line++;
      $rec{line}[$line] = $_;
      $rec{lines}++;
      $found = 1;
      #
      # Query name
      #
      $rec{name} = $1;

      #
      # Query length
      #
      $_ = <$fh>;
      chomp($_);

      $line++;
      $rec{line}[$line] = $_;
      $rec{lines}++;

      while (!(/letters\)/)){
        $_ = <$fh>;
        chomp($_);
        $line++;
        $rec{line}[$line] = $_;
        $rec{lines}++;
      }

      $_=~ s/[a-zA-Z\(\)\s+]//g;

      $rec{Qlen} = $_;
      print ("XX$rec{Qlen}\n");
      $_ = <$fh>;
      chomp($_);
      $line++;
      $rec{line}[$line] = $_;
      $rec{lines}++;
    } # if Query=

   if (/^>/){
     #
     # A  reference hit
     #
     $i++;
     $rec{Qhits}++;
     $rec{Qhit}[$i]{name} = substr($_,1);

     $_ = <$fh>;
     chomp($_);
     $line++;
     $rec{line}[$line] = $_;
     $rec{lines}++;

     while (!(/Length =/)){
       $rec{Qhit}[$i]{name} .= $_;
       $_ = <$fh>;
       chomp($_);
       $line++;
       $rec{line}[$line] = $_;
       $rec{lines}++;
     }


     $rec{Qhit}[$i]{comment} = $rec{Qhit}[$i]{name};
     $rec{Qhit}[$i]{comment} =~ s/\s+/ /g;
     @words = split(/\s+/, $_);
     $rec{Qhit}[$i]{len} = $words[3];
     $rec{Qhit}[$i]{name} =~ s/\s+/ /g;
     @tmp = split(/\s+/,$rec{Qhit}[$i]{name});

     my $len = length($tmp[0]);

     if ($len > $MAXLEN){
       # $tmp = substr($rec{Qhit}[$i]{name}, 0 , $MAXLEN);
       $rec{Qhit}[$i]{name} = substr($tmp[0], 0 , $MAXLEN);
     }
     else{
       $rec{Qhit}[$i]{name} = $tmp[0];
     }
     #
     # Length for a reference hit
     #
     $_=~ s/[a-zA-Z\=\(\)\s+]//g;
     $rec{Qhit}[$i]{len} = $_;
     $_ = <$fh>;
     chomp($_);
       $line++;
       $rec{line}[$line] = $_;
       $rec{lines}++;

    #
    # initialize number of sub-alignments
    #
    $j=0;
   } # if (/^>/){

   if (/Expect/){
     #
     # E-value
     #
     my @words = split(/\s+/, $_);
     $j++;
     $rec{Qhit}[$i]{nseq} = $j;
     $rec{Qali}[$i][$j]{evalue} = $words[8];
     $rec{Qali}[$i][$j]{frame} = "undef";
     $_ = <$fh>;
   }
     print ("XX12 '$rec{Qlen}'\n");

   if (/Identities = (\S+)/){
     print ("XX2 '$rec{Qlen}'\n");
     #
     # Number of identities in alignment
     #
     my @words = split (/\//, $1);
     $rec{Qali}[$i][$j]{mlen} = $words[0];
     $rec{Qali}[$i][$j]{alen} = $words[1];
     if (/Gaps = (\S+)/){
       my @words = split (/\//, $1);
       $rec{Qali}[$i][$j]{gaps} = $words[0];
       # print ("Gaps='$rec{Qali}[$i][$j]{gaps}'\n");
     }

     #
     # percent matched
     #
     $rec{Qali}[$i][$j]{perc_ag} = ($rec{Qali}[$i][$j]{mlen} + $rec{Qali}[$i][$j]{gaps})* 100.0 / $rec{Qali}[$i][$j]{alen};
     $rec{Qali}[$i][$j]{perc_a} = $rec{Qali}[$i][$j]{mlen} * 100.0 / $rec{Qali}[$i][$j]{alen};
     print ("XX$rec{Qlen}\n");
     $rec{Qali}[$i][$j]{perc_q} = $rec{Qali}[$i][$j]{mlen} * 100.0 / $rec{Qlen};

     $_ = <$fh>;
     chomp($_);
     $line++;
     $rec{line}[$line] = $_;
     $rec{lines}++;
   }

   if (/ Frame =/){
     #
     # Read the frame information
     #
     my @words = split(/\s+/, $_);
     $rec{Qali}[$i][$j]{frame} = $words[3];
     #print ("# Frame= '$rec{Qali}[$i][$j]{frame}'\n");
   }

    if (/Strand =/){
      my @words = split(/\s+/, $_);
      if ($words[3] eq "Plus"){
        $rec{Qali}[$i][$j]{frame}="+";
      }
      if ($words[3] eq "Minus"){
        $rec{Qali}[$i][$j]{frame} ="-";
      }

      if ($words[5] eq "Plus"){
        $rec{Qali}[$i][$j]{frame} .="+";
      }
      if ($words[5] eq "Minus"){
        $rec{Qali}[$i][$j]{frame} .="-";
      }
    }

   #
   # Data for the query sequence
   #
   if (/Query: /){
     @words = split(/\s+/, $_);

     if ( ! defined($rec{Qali}[$i][$j]{from})){
       $rec{Qali}[$i][$j]{from} = $words[1];
     }

     $rec{Qali}[$i][$j]{to} = $words[3];     
     
     #
     # Read the Query sequence
     #
     $rec{Qali}[$i][$j]{seq} .= $words[2];
 
     $_ = <$fh>;
     chomp($_);
     $line++;
     $rec{line}[$line] = $_;
     $rec{lines}++;
   } 

   #
   # Data for the Reference sequence
   #
   if (/Sbjct: /){
     #
     # June 20, 2003
     # check for lines like these
     #
     # Sbjct: 1032ELTMFQNDFEKACQAKSEALVLREKSTLERIHKHQEIETKEIYAQRQLLLKDMDLLRGRE 1211

     @words = split(/\s+/, $_);
     if ( ! defined($rec{Dali}[$i][$j]{from})){
       if ($words[1]=~/\D+/){
         print ("her\n");
         @w = split(/\D+/, $words[1]);
	 $words[3] = $words[2];
	 $words[2] = $w[1];
	 $words[1] = $w[0];
       }
       $rec{Dali}[$i][$j]{from} = $words[1];
     }
     if ($words[1]=~/\D+/){
       @w = split(/\D+/, $words[1]);
       $words[3] = $words[2];
       $words[2] = $w[1];
       $words[1] = $w[0];
     }
     $rec{Dali}[$i][$j]{to} = $words[3];

     #
     # Read the reference sequence
     #
     $rec{Dali}[$i][$j]{seq} .= $words[2];
     $_ = <$fh>;
     chomp($_);
     $line++;
     $rec{line}[$line] = $_;
     $rec{lines}++;
   }


} # while ((! eof ($fh)) && !(/^S2:/) ){

  #
  # Put Query and reference sequences into arrays
  #
  for (my $i=1; $i <= $rec{Qhits}; $i++){
    for (my $j=1; $j <= $rec{Qhit}[$i]{nseq}; $j++) {
      @Qseq = ();
      @Dseq = ();
      $Qcnt = 0;

      #
      # Replace X' in Query sequence
      #
#      if ($replace eq "ON"){
#        $Qcnt = $rec{Qali}[$i][$j]{seq} =~ tr/X/X/;
#      }
      @Qseq = split(//, $rec{Qali}[$i][$j]{seq});
      @Dseq = split(//, $rec{Dali}[$i][$j]{seq});
      $k=0;
      while (defined($Qseq[$k])){
        $rec{Qseq}[$i][$j][$k+1] = $Qseq[$k];
        $rec{Dseq}[$i][$j][$k+1] = $Dseq[$k];

        $k++;
      }
#      if ($Qcnt > 0){
        #
        # There were X's in query sequence and they will now be
        # replaced with the original amino/nucleotide letters
        # This can also change the percent identity
        #
#        my $from = $rec{Qali}[$i][$j]{from};

        #
        # Only replace X's for $i == 1 to save time
        #
#        if ($i == 1){
#          repBLAST(\%rec, ,$i, $j, $from, "$idx_file");
#        }
#      }
    }
 }
return (%rec);
} # readBLAST

sub writeBLAST {
  ($rh_rec, $fh) = @_;

  for ($i=1; $i <= $rh_rec->{lines}; $i++){
    print $fh ("$rh_rec->{line}[$i]");
  }
}
#############################################################################
#
# writeBLAST_short a short output format
#
##############################################################################
sub writeBLAST_short{
  ($rh_rec, $fh, $evalue, $mismatch_a, $mismatch_q, $Nhit) = @_;

  if ($Nhit eq "all"){
    $hits1 = $rh_rec->{Qhits};
  }
  else{
    $hits1 = 1;
  }

  if ($rh_rec->{Qhits} == 0){
    print $fh ("$rh_rec->{name}\tnohit\t-\n");
  }
  else{

    for ($i=1; $i <= $hits1; $i++){
      if ($Nhit eq "all"){
        $hits2 = $rh_rec->{Qhit}[$i]{nseq};
      }
      else{
        $hits2 = 1;
      }
      for ($j=1; $j <= $hits2; $j++){
        #
        # evalue criteria ok ?
        #
        if ( ($rh_rec->{Qali}[$i][$j]{evalue} <= $evalue) && 
           ($mismatch_a + $rh_rec->{Qali}[$i][$j]{perc_a} >= 100) &&
           ($mismatch_q + $rh_rec->{Qali}[$i][$j]{perc_q} >= 100) &&
           ($min_mlen <= $rh_rec->{Qali}[$i][$j]{mlen})){
          print  $fh ("$rh_rec->{name}\tHIT\t$rh_rec->{Qhit}[$i]{name}\n");
        }
        else{
          print  $fh ("$rh_rec->{name}\tnohit\t$rh_rec->{Qhit}[$i]{name}\n");
        }
      }
    }
  }
}

#############################################################################
#
# writeBLAST_long a long format
#
##############################################################################
sub writeBLAST_long{
  ($rh_rec, $fh, $evalue, $mismatch_a, $mismatch_q, $min_mlen, $Nhit) = @_;


  if ($Nhit eq "all"){
    $hits1 = $rh_rec->{Qhits};
  }
  else{
    $hits1 = 1;
  }


  if ($rh_rec->{Qhits} == 0){
    print $fh ("$rh_rec->{name}\tnohit\t-\t-\t-\t$rh_rec->{Qlen}\t-\t-\t-\t-\t-\t-\t-\t-\t-\tundef\n");
  }
  else{
    for ($i=1; $i <= $hits1; $i++){

      #
      # Show top hit or all hits
      #
      if ($Nhit eq "all"){
        $hits2 = $rh_rec->{Qhit}[$i]{nseq};
      }
      else{
        $hits2 = 1;
      }

      $short_name = substr($rh_rec->{Qhit}[$i]{name}, 0, 10);
      for ($j=1; $j <= $hits2; $j++){
        #
        # evalue criteria ok ?
        #
        if ( ($rh_rec->{Qali}[$i][$j]{evalue} <= $evalue) &&
           ($mismatch_a + $rh_rec->{Qali}[$i][$j]{perc_a} >= 100) &&
           ($mismatch_q + $rh_rec->{Qali}[$i][$j]{perc_q} >= 100) &&
           ($min_mlen <= $rh_rec->{Qali}[$i][$j]{mlen})){
             print  $fh ("$rh_rec->{name}\tHIT\t$rh_rec->{Qhit}[$i]{name}\t");
        }
        else{
             print  $fh ("$rh_rec->{name}\tnohit\t$rh_rec->{Qhit}[$i]{name}\t");
        }
        print  $fh ("$rh_rec->{Qali}[$i][$j]{mlen}\t");
        print  $fh ("$rh_rec->{Qali}[$i][$j]{alen}\t");
        print  $fh ("$rh_rec->{Qlen}\t");
        printf $fh ("%5.1f\t", $rh_rec->{Qali}[$i][$j]{perc_a});
        printf $fh ("%5.1f\t", $rh_rec->{Qali}[$i][$j]{perc_ag});
        printf $fh ("%5.1f\t", $rh_rec->{Qali}[$i][$j]{perc_q});
        print  $fh ("$rh_rec->{Qali}[$i][$j]{evalue}\t");
        print  $fh ("$rh_rec->{Qali}[$i][$j]{from}\t$rh_rec->{Qali}[$i][$j]{to}\t");
        print  $fh ("$rh_rec->{Qhit}[$i]{len}\t");
        print  $fh ("$rh_rec->{Dali}[$i][$j]{from}\t$rh_rec->{Dali}[$i][$j]{to}\t");
        print  $fh ("$rh_rec->{Qali}[$i][$j]{frame}\n");
      }
    }
  }
}

######################################################################################
#
# replace low complexity regions (X's) with original sequence
#
######################################################################################
sub repBLAST{
  ($rh_rec, $I, $J, $from, $idx_file) = @_;

  my $sequence = "";
  my $repseq = "";
  my @seq = ();   
  my $found = 0;
  $old_perc_a = $rh_rec->{Qali}[$I][$J]{perc_a};
  $old_perc_q = $rh_rec->{Qali}[$I][$J]{perc_q};

  $cmd = "get -n \"$rh_rec->{name}\" -I \"$idx_file\"";
  
  #
  # Read the fasta entry that has been fetched using the program 'get'
  #
  open(TMP, "$cmd |");
  while (! eof (TMP)){
    $line = <TMP>;
    chomp($line);
    if ($line =~ /^>/){
      next;
    }
  $sequence.=$line;
  }

  @seq = split (//, $sequence);
  
  $k=0;
  $skip = 0;
  while (defined($seq[$from + $k - 1 - $skip])){
    $Qaa = $rh_rec->{Qseq}[$I][$J][$from + $k];
    $Daa = $rh_rec->{Dseq}[$I][$J][$from + $k];

    if ($Qaa eq "-"){
      $skip++;
      $k++;
      next;
    }

    if ($Qaa eq "X"){
       my $num = $from + $k - 1 - $skip;
      if ($Daa eq $seq[$num]){
        #
        # Increase the match len  
        #

        $rh_rec->{Qali}[$I][$J]{mlen}++;
        $found++;
      }
    }
    $k++;
  }
  if ($found > 0){
    $rh_rec->{Qali}[$I][$J]{perc_a} = $rh_rec->{Qali}[$I][$J]{mlen} * 100.0 / $rh_rec->{Qali}[$I][$J]{alen};
    $rh_rec->{Qali}[$I][$J]{perc_q} = $rh_rec->{Qali}[$I][$J]{mlen} * 100.0 / $rh_rec->{Qlen};
  }
 return(%rec);
}


sub readBLASTPROF{
  ($fh) = @_;

  my %rec = ();
  my $i = 0;

  while (! eof ($fh)){
    my $j = 2;
    $_ = <$fh>;

    $_=~ s/^\s+//;
    $x=length($_);

    if ($x <=1){next}
    chomp($_);
    @words = split(/\s+/, $_);

    
    if ($words[0]=~/[\d+]/){
        $i++;
        $rec{len}++;

        $rec{resnum}[$i] = $i;
	$rec{seq}[$i] = $words[1];
        $rec{sec}[$i] = ".";
        
        $col = 1;
        while (defined($words[$j])){
          $rec{score}[$i][$col] = $words[$j];
          $col++;
          $j++;
        } 
    }
  }
 return (%rec);
}

sub writeBLAST_comment{
  ($rh_rec, $fh, $evalue, $mismatch_a, $mismatch_q, $min_mlen, $Nhit) = @_;


  if ($Nhit eq "all"){
    $hits1 = $rh_rec->{Qhits};
  }
  else{
    $hits1 = 1;
  }


  if ($rh_rec->{Qhits} == 0){
    print $fh ("$rh_rec->{name}\tnohit\t-\t-\t-\t$rh_rec->{Qlen}\t-\t-\t-\t-\t-\t-\t-\t\n");
  }
  else{
    for ($i=1; $i <= $hits1; $i++){

      #
      # Show top hit or all hits
      #
      if ($Nhit eq "all"){
        $hits2 = $rh_rec->{Qhit}[$i]{nseq};
      }
      else{
        $hits2 = 1;
      }

      $short_name = substr($rh_rec->{Qhit}[$i]{name}, 0, 10);
      for ($j=1; $j <= $hits2; $j++){
        #
        # evalue criteria ok ?
        #
        if ( ($rh_rec->{Qali}[$i][$j]{evalue} <= $evalue) &&
           ($mismatch_a + $rh_rec->{Qali}[$i][$j]{perc_a} >= 100) &&
           ($mismatch_q + $rh_rec->{Qali}[$i][$j]{perc_q} >= 100) &&
           ($min_mlen <= $rh_rec->{Qali}[$i][$j]{mlen})){
            print  $fh ("$rh_rec->{name}\tHIT\t$rh_rec->{Qhit}[$i]{comment}\t");
        }
        else{
          print  $fh ("$rh_rec->{name}\tHIT\t$rh_rec->{Qhit}[$i]{comment}\t");
        }
        print  $fh ("$rh_rec->{Qali}[$i][$j]{mlen}\t");
        print  $fh ("$rh_rec->{Qali}[$i][$j]{alen}\t");
        print  $fh ("$rh_rec->{Qlen}\t");
        printf $fh ("%5.1f\t", $rh_rec->{Qali}[$i][$j]{perc_a});
        printf $fh ("%5.1f\t", $rh_rec->{Qali}[$i][$j]{perc_ag});
        printf $fh ("%5.1f\t", $rh_rec->{Qali}[$i][$j]{perc_q});
        print  $fh ("$rh_rec->{Qali}[$i][$j]{evalue}\t");
        print  $fh ("$rh_rec->{Qali}[$i][$j]{from}\t$rh_rec->{Qali}[$i][$j]{to}\t");
        print  $fh ("$rh_rec->{Qhit}[$i]{len}\t");
        print  $fh ("$rh_rec->{Dali}[$i][$j]{from}\t$rh_rec->{Dali}[$i][$j]{to}\t");
        print  $fh ("$rh_rec->{Qali}[$i][$j]{frame}\n");
      }
    }
  }
}
1;

