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
#use Smart::Comments '###';
require "/home/mizuno_progno/lib/bioinfo.pm";
require "/home/mizuno_progno/lib/general.pm";
require "/home/mizuno_progno/lib/math.pm";

$MAFFT_PATH = '/home/mizuno_progno/freeware/mafft-7.205-with-extensions';
$ENV{ 'MAFFT_BINARIES' } = $MAFFT_PATH . '/binaries';

&find( \&main, @ARGV[ 0 ] );
############################################################
sub main
 {
  my( $path ) = $File::Find::name;

  $path =~ /\.seq$/ || return;
  ### $path
  my $mafft = &basename( $_, '.seq' ) . '.mafft';
  if( ! -f $mafft )
   {
    `$MAFFT_PATH/scripts/mafft --retree 2 --maxiterate 1000 --treeout $_ > $mafft`;
   }
  &parse_mafft( $mafft );
 }
############################################################
sub parse_mafft
 {
  my( $mafft ) = @_;

  my( %RES, $id, $seq );
  open( IN, "<$mafft" );
  while( <IN> )
   {
    chomp;
    if( my( $tmp ) = />(\S+)/ )
     {
      $seq && ( $RES{ $id } = $seq );
      $id = $tmp;
      undef $seq;
     }
    else
     {
      $seq .= $_;
     }
   }
  $seq && ( $RES{ $id } = $seq );
  close IN;

  my $contig = &check_contig( values %RES );
  my $contig_nogap = $contig;
  $contig_nogap =~ s/-//g;
  
  open( OUT, ">$mafft.csv" );
  print OUT join( "\t", 'SEQ_ID', 'SEQ', 'SCORE' ), "\n";
  print OUT join( "\t", $mafft, $contig ), "\n";
  print OUT join( "\t", "$mafft\_nogap", $contig_nogap ), "\n";
  
  foreach my $i ( keys %RES )
   {
    print OUT join( "\t", $i, $RES{ $i }, &nuc_diff( $contig, $RES{ $i } ) ), "\n";
   }
  close OUT;
  &plot_hist( "$mafft.csv" );
 }
############################################################
sub plot_hist
 {
  my( $file ) = shift;

  my $label = &basename( $file );
  open( R, "| R --vanilla" );
  print R <<END_OF_LINE;
dat <- read.table( '$file', sep='\t', header=T, fill=T );
bitmap( '$file.png', 'png256' )
dat <- na.omit( dat )
label <- paste( '$label (n=', length( dat[ , 'SCORE' ] ), ')', sep='' )
hist( dat[ , 'SCORE' ], breaks=seq( 0, 1.0, 0.01 ), main=label, xlab='Distance' );
dev.off()
END_OF_LINE
  close R;
 }
############################################################
sub check_contig
 {
  my( $contig );
  my $len = &max( map{ length( $_ ) } @_ );
  
  for( my $i = 0; $i < $len; $i++ )
   {
    my %CNT;
    map{ $CNT{ substr( $_, $i, 1 ) }++ } @_;
    my @dat = reverse sort{ $CNT{ $a } <=> $CNT{ $b } } keys %CNT;
    $contig .= $CNT{ $dat[ 0 ] } == $CNT{ $dat[ 1 ] } ? '?' : $dat[ 0 ];
   }
  return $contig;
 }
############################################################
sub nuc_diff
 {
  my( $contig, $target ) = @_;

  my( $all, $hit );
  for( my $i = 0; $i < length( $contig ); $i++ )
   {
    my $nuc1 = substr( $contig, $i, 1 );
    my $nuc2 = substr( $target, $i, 1 );
    $nuc1 =~ /[ATGCatgc]/ || next;
    $all++;
    $nuc1 eq $nuc2 || next;
    $hit++;
   }

  return( 1 - $hit / $all );
 }
