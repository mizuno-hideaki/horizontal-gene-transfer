#! /usr/bin/perl

############################################################
#
#  creator: Hideaki Mizuno
#  created: 2010/04/11
#
############################################################
use DBI;
use File::Find;
use File::Basename;
use Smart::Comments '###';
require "/home/mizuno_progno/lib/bioinfo.pm";
require "/home/mizuno_progno/lib/general.pm";
require "/home/mizuno_progno/lib/math.pm";

$CUTADAPT = "/home/mizuno_progno/freeware/cutadapt-1.8.1/bin/cutadapt";
$SEQKIT   = "/home/mizuno_progno/freeware/seqkit";
$CUTOFF    = 70.0;
$READSCORE = $ARGV[ 1 ];
&main( $ARGV[ 0 ] );
############################################################
sub main
 {
  my( $file ) = @_;
  open( LIST, "<$file" );
  open( LOG, ">log.csv" );
  print LOG "CUTOFF\t$CUTOFF\n";
  print LOG "READSCORE\t$READSCORE\n";

  while( <LIST> )
   {
    my @ary = &parse_line( $_ );
	$ary[ 0 ] =~ /\D/ && next;
	$ary[ 1 ] || next;
    my $label = sprintf( "F_%d--R_%d", @ary[ 0, 0 ] );
    my $name  = sprintf( "%d_%s", @ary[ 0, 1 ] );
    $name =~ s/\W/_/g;
    system( "$SEQKIT fq2fa $label.fastq > $label.fasta" );
    system( "/home/mizuno_progno/freeware/ncbi-blast-2.2.29+-src/c++/ReleaseMT/bin/makeblastdb -in $label.fasta -dbtype nucl" );
    system( "/home/mizuno_progno/freeware/ncbi-blast-2.2.29+-src/c++/ReleaseMT/bin/blastn -query ../BovB_amplified_Region.fas.txt -outfmt 6 -db $label.fasta -out $label.res -word_size 4 -max_hsps 1 -num_alignments 100000" );
    system( sprintf( "$CUTADAPT -g %s $label.fasta      -e 0.2 --overlap 10 --discard-untrimmed > $label.tmp1.fasta", $ary[ 3 ] ) );
    system( sprintf( "$CUTADAPT -a %s $label.tmp1.fasta -e 0.2 --overlap 10 --discard-untrimmed > $label.tmp2.fasta", &complementary_seq( $ary[ 5 ] ) ) );
    system( sprintf( "$CUTADAPT -g %s $label.fasta      -e 0.2 --overlap 10 --discard-untrimmed > $label.tmp3.fasta", $ary[ 5 ] ) );
    system( sprintf( "$CUTADAPT -a %s $label.tmp3.fasta -e 0.2 --overlap 10 --discard-untrimmed > $label.tmp4.fasta", &complementary_seq( $ary[ 3 ] ) ) );
    &load_blast_res( "$label.res" );
    print LOG join( "\t", $label, $name, '+', &output_seq( $label, $name, '+' ),
                                         '-', &output_seq( $label, $name, '-' ) ), "\n";
    unlink( "$label.tmp1.fasta", "$label.tmp2.fasta", "$label.tmp3.fasta", "$label.tmp4.fasta" );
   }
  close LOG;
  close LIST;
 }
############################################################
sub output_seq
 {
  my( $label, $name, $strand ) = @_;
  my( @log );
  open( OUT, ">>$name.csv" );
  open( SEQ, ">>$name.seq" );
  
  $strand eq '+' && print OUT join( "\t", 'SEQ_ID', 'QV', 'BLAST_SCORE', 'STRAND', 'SEQUENCE' ), "\n";
  open( IN, sprintf( "<$label.tmp%d.fasta", $strand eq '+' ? 2 : 4 ) );

  while( ! eof IN )
   {
    my $header = <IN>;
    my $seq    = <IN>;
	$log[ 0 ]++;
    chomp( $seq );
    $strand eq '-' && ( $seq = &complementary_seq( $seq ) );
    my @ary    = &parse_line( $header );
    my @label  = $ary[ 0 ] =~ /^>(.+)\s(.+)\s/;
    $label[ 1 ] < $READSCORE && next;
    $log[ 1 ]++;
    $BLAST{ $label[ 0 ] }{ 'SCORE' } < $CUTOFF && next;
    $log[ 2 ]++;
    $BLAST{ $label[ 0 ] }{ 'LENGTH' } < 1000 && next;
    $log[ 3 ]++;
    $BLAST{ $label[ 0 ] }{ 'STRAND' } eq $strand || next;
    $log[ 4 ]++;
    length( $seq ) < 1300 && next;
    $log[ 5 ]++;
    length( $seq ) > 3500 && next;
    $log[ 6 ]++;
    print OUT join( "\t", @label, $BLAST{ $label[ 0 ] }{ 'SCORE' }, $strand, $seq ), "\n";
    print SEQ '>', $label[ 0 ], "\n";
    print SEQ $seq, "\n";
    $DAT{ $label[ 0 ] }++;
    $DAT{ $label[ 0 ] } > 1 && die;
   }
  close IN;

  close SEQ;
  close OUT;
  
  return @log[ 0 .. 6 ];
 }
############################################################
sub load_blast_res
 {
  my( $blast ) = @_;
  ### $blast

  open( BLAST, "<$blast" );
  $_ = <BLAST>;
  while( <BLAST> )
   {### %
    my @ary = &parse_line( $_ );
    $BLAST{ $ary[ 1 ] }{ 'SCORE' }  = $ary[ 2 ];
	$BLAST{ $ary[ 1 ] }{ 'LENGTH' } = $ary[ 3 ];
    $BLAST{ $ary[ 1 ] }{ 'STRAND' } = $ary[ 8 ] < $ary[ 9 ] ? '+' : '-';
   }
  close BLAST;
 }

