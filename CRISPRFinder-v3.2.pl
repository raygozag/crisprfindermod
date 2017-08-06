#!/usr/bin/env perl
# crispr finder version 3
# Ibtissem GRISSA 24/05/2006
# this script allows to find crisprs and to store the result in the files
# Modified by Christine Drevet to work with EMBOSS Version 5 (01/12/2009)
# Modified by Luis M. Rodriguez-R <lrr@cpan.org> to improve portability and code readability (Feb 2011 and Jun 2011)
# Modified by Lee S. Katz to produce gff files (Feb-March 2011)

# CPAN modules 
use strict;
use Class::Struct;   #charge le module qui construit les struct
use warnings;
use Bio::SeqIO;
use File::Copy;
use File::Basename;
use Data::Dumper;
use Cwd;
# TODO add URLescape module in case the GFF attributes get messy
#local modules

# set parameters for the Vmatch program
#$ENV{'LD_LIBRARY_PATH'} = '.';

if(@ARGV<1)
  {
    printhelpbasic($0);
    exit 1;
  }

  if ($ARGV[0] eq '-help')
  {
    printhelpall($0);
    exit 0;
  }

  if ($ARGV[0] eq '-v')
  {
    printversion($0);
    exit 0;
  }


my $userfile = $ARGV[0];
if (!$ARGV[1]) { $ARGV[1] = './';}
my $ResultDir = $ARGV[1];

my $seqIO = Bio::SeqIO->new(-format=>'Fasta', -file=>$userfile);
my $inputfileCount=0;
my ($seq,$inputfile);
my $basename = basename($userfile);
my @outdir = split /\./, $basename;
my $outdir = $outdir[0];
my @set = ('0' ..'9', 'A' .. 'F');
my $strk = join '' => map $set[rand @set], 1 .. 8;
# LK
#if(-d $ResultDir){
#  die "Output directory $ResultDir already exists. Remove the directory before continuing.\nSuggestion:\nrm -rf $ResultDir/\n";
#}
mkdir $outdir unless -d $outdir;

  struct Rep => {
    Pos1  => '$',
    Length  => '$',
    DRseq => '$',
  };

my $totalNumberOfCrisprs=0; #LK
while($seq = $seqIO->next_seq){
  $inputfileCount++;
  $inputfile = $basename;
  #print getcwd()."\n";
  my $seqO = Bio::SeqIO->new(-file=>">$inputfile", -format=>'Fasta');
  $seqO->write_seq($seq);
  
  my $nbrcris = 0;
  
  # execute the mkvtree function to indexate the input sequence file
  my $indexname = $inputfile;
  my @indexname = split (/\./,$indexname);
  #print STDERR "Sequence number $inputfileCount..";
  #print STDERR "\n*** your results files will be in the $outdir-$inputfileCount directory ***\n"; 
  
  callmkvtree($inputfile,$indexname);  
  # execute the first vmatch search
  # this line
  my @vmatchoptions = qw(-l 23 25 100000 -e 1 -s leftseq -evalue 1 -absolute -nodist -noevalue -noscore -noidentity -sort ia -best 1000000 -selfun sel392.so 100);
  
  push(@vmatchoptions,$indexname);
  push(@vmatchoptions, " > /tmp/vmatch_result${strk}.txt");
  # makesystemcall("./vmatch " . join(' ',@vmatchoptions));
  
  makesystemcall("vmatch " . join(' ',@vmatchoptions)); #LK
  makesystemcall("python2.7 /home/qiime/Downloads/process.py /tmp/vmatch_result${strk}.txt /tmp/vmatch_result${strk}2.txt");
  my @rep = trans_data("/tmp/vmatch_result${strk}2.txt");
  # fill in tabseq table : table containing the begin and end positions 
  # of sequences candidates ton contain a crispr
  my @temp = split(/\./, $inputfile);
  my $RefSeq = $temp[0];
  #my $OneSpacerCris_nbr=0;
  if($#rep >=0)
  {
     my($nbrcris, $OneSpacerCris_nbr)=  write_clusters($RefSeq,@rep);
     create_recap($RefSeq, $nbrcris, $OneSpacerCris_nbr,$ResultDir);
     $totalNumberOfCrisprs+=$OneSpacerCris_nbr; #LK
      `rm -Rf $outdir-$inputfileCount &>/dev/null`;
  `mv $ResultDir/$outdir $ResultDir/$outdir-$inputfileCount 2>&1 >/dev/null`; #LK     
     if($nbrcris>0 || $OneSpacerCris_nbr>0){
  
     #`cp ../$inputfile $ResultDir/$outdir-$inputfileCount/$inputfileCount.fasta `;
     }
     
     rmdir "$outdir-$inputfileCount"; # Only if empty
  }
  else{ create_recap($RefSeq, 0, 0,$ResultDir); }
  $totalNumberOfCrisprs+=$nbrcris; #LK

  
  #unlink <..\/crispr_result_*>;
  #unlink <..\/seq_v*>;
  #unlink <..\/$inputfile.*> ;
  #unlink <..\/spacersseq_v*>;
  #unlink '../stdout';
  #unlink '../vmatch_result.txt';
  #unlink '../alignDR_Spacer.needle';
  #unlink '../DR';
  #unlink <*_CRISPRs> ;
  #unlink "../$inputfile";
  
  
  #rmdir $outdir;
  
} # end while (my $seq = $seqIO->next_seq)
print STDERR "\n"; # LK

# LK
my $gffFilename=makeGff($ResultDir,$inputfile);
my $seqIdPrefix=(split(/\//,$userfile))[-2]; # prefix is just before assembly.fasta
#print getcwd()."\n";
my $fastaFilename=makeFasta($gffFilename,$userfile,$seqIdPrefix);
print STDERR "CRISPRfinder found $totalNumberOfCrisprs CRISPRs\n";
# // LK

#TODO make a GBROWSE file too
# http://nbase.biology.gatech.edu/gbrowse2/general_help.html#upload

#------------------------------------------------------------------------------
# /CRISPRFinder functions

# LK
sub makeFasta{
  my ($gffFilename,$genomeFilename,$seqIdPrefix)=@_;
  my($filename,$dir)=fileparse($gffFilename);
  my $basename=basename($filename,'.gff');
  my $fastaFilename="$dir/$basename.fasta";

  # load contigs/chromosomes into a hash
  my %seq;
  my $seqin=Bio::SeqIO->new(-file=>$genomeFilename,-format=>"fasta");
  while(my $seq=$seqin->next_seq){
    $seq{$seq->id}=$seq;
  }

  open GFF, $gffFilename or die("Could not open the GFF file because $!\n");
  my $seqout=Bio::SeqIO->new(-format=>"fasta",-file=>">$fastaFilename");
  my $crisprCounter=0;
  while(<GFF>){
    s/^\s+|\s+$//g; # trim whitespace
    next if(/^#/ || /^$/); # skip comments or blank lines

    my($seqid,$source,$type,$start,$stop,$phase,$strand,$score,$attributes)=split(/\t/);
    next if($type ne 'CRISPR'); # skip non-crispr lines
    
    # write the sequence
    my $subseq=$seq{$seqid}->subseq($start,$stop);
    my $id=join("_",$seqIdPrefix,$seq{$seqid}->id,++$crisprCounter);
      $id=~s/[\|\.]//g; # remove some characters from the ID
    my $desc="Length=".($stop-$start+1)." Start=$start Stop=$stop Contig=".$seq{$seqid}->id."";
    my $seqObj=Bio::Seq->new(-id=>$id,-seq=>$subseq,-desc=>$desc);
    $seqout->write_seq($seqObj);
  }
  close GFF;

  return $fastaFilename;
}
#LK
# parse result files and create a GFF file for all
sub makeGff{
  my($ResultDir,$inputfile)=@_;
  my @fastaExtensions=qw(.fasta .fna .faa .mfa .fa .fas .ffn);
  my $inputBasename=basename($inputfile,@fastaExtensions);
  my @dir=<$ResultDir/$inputBasename*>;

  # parse
  my $GFFstr="";
  open GFF,">$ResultDir/$inputBasename.gff" or die "Could not open output GFF in $ResultDir/ because $!\nI will not create the GFF file.\n";
  for my $dir(@dir){
    my @spacerFile=<$dir/Spacers_*>;
    for my $spacerFile(@spacerFile){
      my($spacerNumber);
      if($spacerFile=~/_(\d+)$/){
        $spacerNumber=$1 ;
      } else {
        warn "Warning: Cannot understand result file $spacerFile. Skipping.\n";
        next;
      }

      # transform the report file to a GFF line
      my $reportFile="$dir/$inputBasename"."_PossibleCrispr_$spacerNumber";
      $reportFile="$dir/$inputBasename"."_Crispr_$spacerNumber" unless -e $reportFile;
      unless(-e $reportFile){
         warn "Impossible to find details for spacer $spacerNumber";
	 next;
      }
      my $feature=reportToGff($reportFile);

      # add onto the GFF
      $GFFstr.="$feature\n" if defined $feature;
    }
  }
  $GFFstr="##gff-version 3\n$GFFstr"; # add in the definition line
  print GFF $GFFstr;
  close GFF;
  #print STDERR "GFF results are  in $ResultDir/$inputBasename.gff\n";

  return "$ResultDir/$inputBasename.gff";
}
#LK
# transform a given report file to a line for a GFF file
# returns a string TODO a feature object
# http://www.bioperl.org/wiki/GFF3
# http://doc.bioperl.org/releases/bioperl-current/bioperl-live/Bio/SeqFeature/Generic.html
sub reportToGff{
  my($reportFile)=@_;
  #GFF columns
  my($seqid,$start,$end,$strand,$attributes,%tag);
  my($source,$type,$score,$phase)=("CRISPRfinder","CRISPR",".",".");
  $type="PossibleCRISPR" if $reportFile =~ /PossibleCrispr_\d+$/;
  
  my($featureId,$gffLine,@spacerInfo,$DRsequence);
  return unless -e $reportFile;
  # get data from the report file
  my $onSpacersList = 0;
  open(REPORT,"<",$reportFile) or die("Could not open the report file $reportFile because $!\n");
  while(<REPORT>){
    my($spacerStart,$spacerEnd,$spacerSequence);
    if(/Id: (\S+)/){
      $seqid=$1;
    }
    elsif (/Crispr.+?(\d+)\s+Crispr\D+(\d+)/){
      $start=$1;
      $end=$2;
      my $nextLine=<REPORT>;
      $nextLine=~s/://g; # remove colons
      my @field=split(/\s+/,$nextLine);
      shift(@field); # remove the hash sign
      for(my $i=0;$i<@field;$i+=2){
        $attributes.="$field[$i]=$field[$i+1];";
        $tag{$field[$i]}=$field[$i+1];
      }
    }
    elsif (/Spacer_begin_position/){
      $onSpacersList = 1;
    }
    elsif ($onSpacersList){
      next if /^#/;
      next if /^\s*$/;
      my @field=split(/\s+/,$_);
      shift(@field);
      my $name=join("_","spacer",$field[0],$field[1]);
      push(@spacerInfo,{
        start=>$field[0],end=>($field[0]+$field[1]-1),
        type=>"CRISPRspacer",
        attributes=>"sequence=$field[2];name=$name",
        }
      );
    }
  }
  close REPORT;

  # derivations of the data
  if($start<$end){
    $strand='1';
  } else{
    $strand='0';
  }
  $featureId=join("_",$seqid,$start,$end);
  $attributes.="name=$featureId;ID=$featureId";
  $DRsequence=$tag{DR}; delete($tag{DR});

  # Construct the gff line with 1) the entire CRISPR and then 2)DRs interspersed with spacers
  $gffLine=join("\t",$seqid,$source,$type,$start,$end,$score,$strand,$phase,$attributes)."\n";
  my $DRstart=$start;
  for my $s (@spacerInfo){
    # DRs surround spacers.  
    $gffLine.=join("\t",$seqid,$source,"CRISPRdr",$DRstart,($DRstart+$tag{DR_length}-1),$score,$strand,$phase,"Parent=$featureId")."\n";
    # spacer
    $$s{attributes}.=";Parent=$featureId";
    $gffLine.=join("\t",$seqid,$source,$$s{type},$$s{start},$$s{end},$score,$strand,$phase,$$s{attributes})."\n";

    # update DR for the next iteration
    $DRstart=$$s{end}+1;
  }
  # one more DR
  $gffLine.=join("\t",$seqid,$source,"CRISPRdr",$DRstart,($DRstart+$tag{DR_length}-1),$score,$strand,$phase,"Parent=$featureId")."\n";

  return $gffLine;
}

sub compare
{
  my ($val1, $val2) = @_;
  return ($val1<=$val2*0.5) || ($val1>=$val2 * 1.5) ? 1 : 0;
  #if( ($val1<=$val2*0.5) || ($val1>=$val2 * 1.5) ){return 1;} else {return 0;}
}

sub extractcrispr
{
  my($seqfile,@spacers) = @_;
  my($count,$crisprfile,$seq,$spacerseq, $beg, $end, $id);

	$crisprfile = "spacers".$seqfile;
	open WRITER,"> $crisprfile" or die "The file cannot be edited !\n";
	my $in = new Bio::SeqIO(-format => 'Fasta',-file => $seqfile);
	$seq = $in->next_seq;

	$count=0;
	while($count<= $#spacers)
	{
		$beg = $spacers[$count]->Pos1;
		$end = $beg + $spacers[$count]->Length-1;
		$spacerseq = $seq->subseq($beg, $end);
		my $u = $count+1;
		$id = "spacer".$u;
		print WRITER ">".$id."\n";
		print WRITER $spacerseq."\n";
		$count++;
	}
	close WRITER;
}
#------------------------------------------------------------------------------
sub definespacers
{  #function that returns a structure of hypothetical crisprs definition (the structure spacers )
  my($DRlength, $crisprfile,$indexname)= @_;
  my(@lines, $count, @temp,$posdeb, $posend, $len,$i,$refFalsSpacers, @spacers_correct);
  my @spacers = ();
  my $nbspacers = 0;

my $mism = 0;
my $Toterr = 0;
my $simDRs = 0; # indicates if all DRs are different
  open(FD, $crisprfile)  or die "Error in opening the input file : $crisprfile!";
  $i = 0;	#counter for the number of spacers we find
  my $repCount = 0;
  my $repCurrent = 0;
  while(my $line = <FD>)
  {
     if($line=~/^# HitCount:\s+(\d+)/i){
        $repCount = $1+0;
     }
     next if $line =~ /^#/;
     $line = trim($line);
     next unless $line;
     next if $line =~ /^Start/;
     $repCurrent++;
     @temp = split(/ +/, $line);
     # store the first DR start end position
     unless(defined $posend){
        $posend = $temp[1];
	next;
     }
     $posdeb = $temp[0];
     $len = $posdeb - $posend -1 ;
     if($len < ($DRlength*0.6) )
     {
	$mism = $temp[$#temp-1] eq "." ? 0 : $temp[$#temp-1]/$DRlength;
     }
     else
     {
	if( $repCurrent == $repCount )
	{
		if($len > ($DRlength*2.5))
		{
			$posend = $temp[1];
		} 
		else
		{
			@spacers = add_spacer($posend,$len,$i,@spacers);
			$simDRs = 1 if $temp[$#temp-1] eq '.';
			$Toterr = $Toterr + $mism;
			$posend = $temp[1];
			$i++;
		}
		$mism = $temp[$#temp-1] eq "." ? 0 : $temp[$#temp-1]/$DRlength;
	}
	else
	{
		@spacers = add_spacer($posend,$len,$i,@spacers);
		$simDRs = 1 if $temp[$#temp-1] eq '.';
		$Toterr = $Toterr + $mism;
		$posend = $temp[1];
		$mism = $temp[$#temp-1] eq "." ? 0 : $temp[$#temp-1]/$DRlength;
		$i++;
	}
     }
   }
   close FD;
   
   $Toterr = ($Toterr + $mism)/($i+1);
   if($Toterr > 0.2)
   {
   	@spacers = ();
	my @false_spacers = ();
	$refFalsSpacers = \@false_spacers;
   }
   
   if($#spacers >=0)
   {
	extractcrispr($indexname,@spacers);
	($refFalsSpacers, @spacers_correct) = check_cris_div($DRlength,$indexname, @spacers);
   }
   else 
   {
   	my @false_spacers = ();
   	$refFalsSpacers = \@false_spacers;
   }
   
   return($simDRs,$refFalsSpacers, @spacers_correct);
}

#------------------------------------------------------------------------------
# this function checks the spacers length (comparison with 2.5 DR)
# deletes the wrong spacers in the begining 
# returns an array of wrong spacers begin positions 
sub check_cris_div
{
  my($DRlength,$indexname, @spacers) = @_;
  my @false_spacers = ();
  my($i,$j,$comp,$begtest,$spLen,$spPos, @spacers_new);
  @spacers_new = ();
  $j = 0;
  $begtest = 1;
  $comp = $DRlength * 2.5;
	my $reffal_spacers;

  for($i=0; $i<= $#spacers; $i++)
  {
		$spLen = $spacers[$i]-> Length;
		if( ($spLen > $comp) )
		{
	  		if($i == $#spacers)
	  		{
				if($#spacers > 0)
				{
					my $crisprfile = "spacers".$indexname;
					open(FD, $crisprfile)  or die "Error in opening the input file : $crisprfile!";
					my @lines=<FD>;
					close(FD);
					open WRITER,"> TmpSpacers" or die "The file TmpSpacers cannot be edited !\n";
					my $k = 0;
					while($k<$#lines -1)
					{
		  				print WRITER $lines[$k];
		  				$k++;
					}
					close WRITER;
					rename("TmpSpacers",$crisprfile);
					# not good : @spacers_new is not the solution!!
					($reffal_spacers,@spacers_new)=check_cris_div($DRlength,$indexname, @spacers_new);
				}
	  		}
	  		else
	  		{
	    		if(!$begtest)
	    		{
					push(@false_spacers,$spacers[$i]-> Pos1);
					$spPos = $spacers[$i]-> Pos1 - 1;
					@spacers_new = add_spacer($spPos,$spLen,$j,@spacers_new);
					$j++; $begtest = 0;
				}
	    		else
	    		{
					my $crisprfile = "spacers".$indexname;
					open(FD, $crisprfile)  or die "Error in opening the input file : $crisprfile!";
					my @lines=<FD>;
					close(FD);
					open WRITER,"> TmpSpacers" or die "The file TmpSpacers cannot be edited !\n";
					my $k = 2;
					while($k<=$#lines)
					{
		  				print WRITER $lines[$k];
		  				$k++;
					}
					close WRITER;
					rename("TmpSpacers",$crisprfile);
	    		}
			}
		}
		elsif($spLen <= $comp)
		{
			$spPos = $spacers[$i]-> Pos1 - 1;
			@spacers_new = add_spacer($spPos,$spLen,$j,@spacers_new);
			$j++; $begtest = 0;
	   }
   }
	if(!$reffal_spacers){$reffal_spacers = \@false_spacers;}
	return($reffal_spacers,@spacers_new);
}
#------------------------------------------------------------------------------
sub add_spacer
{
  my($pos,$len,$i,@spacers)=@_;
  $pos = $pos+1;
  $spacers[$i] = Rep->new();
  $spacers[$i]-> Pos1($pos);
  $spacers[$i]-> Length($len);

  return @spacers;
}
#------------------------------------------------------------------------------
# created 26/10/2006
# use of clustalw in alignement of multiple spacers
# modif 16/02/2007
# eval to go on if clustalw errors
sub checkspacersAlign
{
 my($file) = @_;
 my $ident;
 eval
 {
 	#$ENV{CLUSTALDIR} = '/home/username/clustalw1.8/';
 	use Bio::Tools::Run::Alignment::Clustalw;
	#open STDOUT, "|tee stdout >/dev/null 2>&1";
	my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
	my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);

	#  Pass the factory a list of sequences to be aligned.	
   my $aln = $factory->align($file); # $aln is a SimpleAlign object.
  	$ident = $aln->percentage_identity;
 };
 
 if($@)
 {
    # Une erreur s'est produite...
	$ident = 100;
 }

	if($ident >= 60) {return 0;} else {return 1;}
}

#------------------------------------------------------------------------------
sub checkspacers
{
  my($len,$file) = @_;
  my($check,@temp);
  # execute the mkvtree function to indexate the input sequence file
  callmkvtree($file,$file);
  
  my $err = int($len /4);
 
  my @vmatchoptions = ("-l", $len, "-e $err", "-nodist","-noevalue", "-noscore", "-noidentity");
  push(@vmatchoptions,$file);
  push(@vmatchoptions, " > check_spacer");
    
  # makesystemcall("./vmatch " . join(' ',@vmatchoptions));
  makesystemcall("vmatch " . join(' ',@vmatchoptions)); #LK
  $check = `wc -l check_spacer`;
  @temp = split(/ +/, $check);
  
   if($temp[0]>=2)
   {
		my $rep = $temp[0];
      $check = `wc -l $file`;	# count the number of spacers *2
  		@temp = split(/ +/, $check);
		my $spac = $temp[0];
		if($rep >= $spac/4) { return 0; }else {return 1;}
   } 
   else {return 1;}
}
#------------------------------------------------------------------------------
sub extractsequence
{
  # extract a genomic sequence from a well defined position till a predefined length
  #extractsequence(FILENAME, startPOS1, endPOS1, startPOS2, endPOS2,...)
  # the output is files containig each hypothetical crispr loci
  my(@seqencesdefinition) = @_;
  my($count,$seqfile, $number,@subselectoptions); 
  for($count = 0; $count < $#seqencesdefinition; $count+=2)
  {
		# extending the sequences delimiters and testing that they doesn't overflow
		#positions cases
    	if($seqencesdefinition[$count] - 500 > 0){$seqencesdefinition[$count]=$seqencesdefinition[$count] - 500;} 
    	else{$seqencesdefinition[$count]=0;}
		if($seqencesdefinition[$count+1]+500> $seq->length()){$seqencesdefinition[$count+1]=$seq->length()-1;}
		else{ $seqencesdefinition[$count+1] = $seqencesdefinition[$count+1] + 500;}
    	$number = ($count+2)/2;
    	$seqfile = "seq_v" . $number;
    	@subselectoptions = ("-range", $seqencesdefinition[$count], $seqencesdefinition[$count+1], $inputfile, ">", $seqfile);
    	# makesystemcall("./vsubseqselect " . join(' ',@subselectoptions));
    	makesystemcall("vsubseqselect " . join(' ',@subselectoptions)); #LK
   }
   return @seqencesdefinition;
}

#------------------------------------------------------------------------------
# transform output file to a structure
sub trans_data
{
  my $file = $_[0];
  my $elem_nbr = 0;
  my @repetitions = ();
  my @temp;
  my $line;
  my $ln=0;
 open(FD, $file)  or die "Error in opening the input file : $!";
 while( <FD> ) 
 {
     next if /^(\s)*$/;  # skip blank lines
     chomp;              # remove trailing newline characters
		$line = $_;
		$line =~ s/>/> /g;
      @temp = split(/ +/, $line);
      $ln=$ln+1;
      if($temp[0] =~ /^\>/) 
      {
			$repetitions[$elem_nbr] = Rep->new();
			$repetitions[$elem_nbr]-> Length($temp[1]);
			$repetitions[$elem_nbr]-> Pos1($temp[2]);
     	} 
     else
     {
			if(defined $temp[0]){
			$repetitions[$elem_nbr]-> DRseq($temp[0]."aaaaaaaaaa");
			$elem_nbr++;
			}
     }
  }
close(FD);
return @repetitions;
}

#------------------------------------------------------------------------------
# find clusters and check DRs ---- 
sub write_clusters{

  my ($RefSeq,@rep) = @_;
  my @tabseq = find_clusters(@rep);

  my($count, $seqbeg, $seqend,$i,$j, $DR, $DRocc, $nbrcris,$DRlength);
  $i=0;
  $j=0;
  $nbrcris = 0;
  my $OneSpacerCris_nbr = 0;
 my $modf=-1;
 my @modifcode =();
  # extract all the candidate sequences and modify the seq delimiters
  @tabseq = extractsequence(@tabseq);	

  # --- analyze the candidate DR in each extracted sequences
  for($count=1; $count<= $#tabseq/2+1; $count++)
  {
		my @DR_cand = ();	# store all the DR candidates
		my @DR_cand_occ = ();	# store the possible DRs after comparing the number of occurences
		$seqbeg = $tabseq[2*$count-2]; #start position of the locus file
		$seqend = $tabseq[2*$count-1]; #end position of the locus file

		# -- store the length of the first DR candidat
		# --------------------
		$DR = $rep[$i]->DRseq;
		push(@DR_cand, $DR);

		my $occ_DR_max = DR_occ_rev($DR, @rep);
		push(@DR_cand_occ,$DR);
		$i++;
		# --------------------

		while(($i<= $#rep) && ($rep[$i]->Pos1 <= $seqend) )
		{
			$DR = $rep[$i]->DRseq;
			if(!check_DR($DR, @DR_cand))
			{
				push(@DR_cand, $DR);
				$DRocc = DR_occ_rev($DR, @rep);  #number of occurences (rev comp included)
				# find the maximal number of occurences (if equal, take the minimal length)
				if($occ_DR_max < $DRocc)
				{
					$occ_DR_max = $DRocc; 
					@DR_cand_occ=(); 
					push(@DR_cand_occ,$DR);
				}
				else
				{
			  		if($occ_DR_max == $DRocc)
			  		{
						if( $DR ne $DR_cand_occ[0] ){push(@DR_cand_occ,$DR);}

			  		}
				}
			}
			$i++;
		}

   	my $indexname = "seq_v".$count;
   	my $crisprfile = "crispr_result_".$count;   
   	my $bestscore = 100000; 
   	my $bestindex = 0;
   	my $c=0;  # counter of the DRs consensus
   	my $TotDRerr=0;my $TotDRerr_best=0;
   	while($c <= $#DR_cand_occ)
   	{
			my $score;
			# case we have more than one consensus for the DR
			if($#DR_cand_occ > 0)
			{
				# we have c DRfile and c crispr_result files
				$crisprfile = "crispr_result_".$count."_".$c;
			}	
			$DRlength = length($DR_cand_occ[$c]);
   		my $err = int($DRlength/3);
   		my @fuzznucoptions = ("-sequence", $indexname, "-pattern", $DR_cand_occ[$c], "-pmismatch", $err, "-outfile",$crisprfile);
		push @fuzznucoptions, "-auto";
   		makesystemcall("fuzznuc " . join(' ',@fuzznucoptions));
		if($#DR_cand_occ > 0)
		{
			#compare the crispr analysis based on the possible DRs
			($score,$TotDRerr) = compute_Crispr_score($crisprfile,length($DR_cand_occ[$c]));
			if($score <= $bestscore){$bestscore = $score; $bestindex = $c;$TotDRerr_best=$TotDRerr;}
		}
		else
		{
			($bestscore,$TotDRerr_best) = compute_Crispr_score($crisprfile,length($DR_cand_occ[0]));
		}
		$c++;
   	} # boucle while
   	# take only the good DR seq and crispr_result seq (rename it to be good!)
   	if($#DR_cand_occ > 0)
   	{
   		$crisprfile = "crispr_result_".$count;
			rename($crisprfile."_$bestindex", $crisprfile);
			$DR = $DR_cand_occ[$bestindex];
   	}
   	else
	{
		$DR = $DR_cand_occ[0];
	}

   	my ($crisOK,$simDRs,$CrisprBeg, $CrisprEnd,$RefSpacersH,$nbspacers,$RefFalsSpacers) =
			Find_theCrispr($indexname,$DR,$count,$seqbeg,$seqend,$crisprfile);
	if($crisOK)
	{
		if($TotDRerr_best > 0.2){$crisOK = 0;}
	} 
   
   	# compute the number of crisprs in the whole sequence
   	$nbrcris = $nbrcris + $crisOK;
	
   	# if the current subsequence contains a crispr, it is stored
   	if($crisOK)
   	{
			my $Crispr_file = fill_in_crisprfile($simDRs,$RefSeq,$ResultDir, $nbrcris, $CrisprBeg, $CrisprEnd,$DR, $nbspacers, "spacers".$indexname,$RefFalsSpacers, %$RefSpacersH);
			my @FalsSpacers = @$RefFalsSpacers;
			if($#FalsSpacers != -1){$modf=$#FalsSpacers;}
			if($nbspacers <= 1 || $simDRs == 0 ){$OneSpacerCris_nbr ++;}
			if($modf != -1)
			{
		 		chdir($ResultDir."/".$RefSeq);
	    		my $hyp_cris_nbr = $OneSpacerCris_nbr;
		 		my $crisprs_nbr = $nbrcris;
		 		#--------------------------------------
		 		while($modf >=0)
		 		{
	  		 		if($crisprs_nbr - $hyp_cris_nbr >= 1)
	  		 		{
			  			my $s = 1;
			 		 	while($s<=$crisprs_nbr)
			 		 	{
			     			$inputfile = "$RefSeq"."_Crispr_".$s;
			     			if(-e $inputfile)
			     			{
							($s,$crisprs_nbr,$OneSpacerCris_nbr) =
								modify_files($inputfile,$ResultDir,
									$RefSeq,$s,$crisprs_nbr,
									$OneSpacerCris_nbr,$modf);
								$nbrcris = $crisprs_nbr;
			     			}
			     			$s++;
			 		 	}
			 		}
	 	   		$modf--;
				}
			}
			$modf = -1;

   	}
 } #boucle for

return ($nbrcris, $OneSpacerCris_nbr);
}

#------------------------------------------------------------------------------
sub modify_files
{
  my($inputfile,$ResultDir,$RefSeq,$rank,$crisprs_nbr,$HYPCris_nbr,$modf)=@_;
  my(@lines,@temp,$nb_spacers,@nb_div, $newCrisNbr, @divi);

  my $dir = $RefSeq."/";
  chdir($dir);

  open(FD, $inputfile)  or die "Error in opening the file:$inputfile!";
  @lines=<FD>;
  close(FD);

  $newCrisNbr = $crisprs_nbr;
  @temp = split(/:/, $lines[15]);
  my @temp2 = split(/ +/, $temp[1]);
  my $begPos = $temp2[1];
  my $endPos = $temp[2];
  @temp = split(/:/, $lines[16]);
  @temp2 = split(/ +/, $temp[2]);
  my $DRlength = $temp2[1];
  $nb_spacers = $temp[3];
  my $j = 20+2*$nb_spacers;

  if($lines[$j] =~ '#'){}
  else{
	# case we have divisions
	@temp = split(/:/, $lines[$j]);
	@nb_div = split(/ +/, $temp[1]);
	if($modf != -1) { @divi = @temp; shift @divi;shift @divi;pop @divi; }
		$newCrisNbr += scalar(@nb_div);
# hypothese : there is only one division per CRISPR
	my $line = find_Sp_lines($j,$inputfile,$nb_div[0]);
	my($file_new, $file_old,$l1,$l2,$id,$CrisprBeg,$CrisprEnd,$nbspacers,$Spfile_old,$Spfile_new);
	$file_old = $inputfile;
	$Spfile_old = "Spacers_".$rank;

my ($nbsp1,$nbsp2);
# ----------file1-----------------------
$l1 = 20;
$l2 = $line ;
#print "rank = $rank inputfile : $inputfile\n-------------\n";
$id = $rank;
$Spfile_new = "Spacers_test_".$id;
$CrisprBeg = $begPos;
$CrisprEnd = $nb_div[0];
$CrisprEnd = $CrisprEnd."\n";
$nbsp1 = ($l2 - $l1)/2;
if($nbsp1>=1){
if($nbsp1 >= 3){$file_new = "$inputfile"."_test_".$id;}
else{$file_new = "$inputfile"."_test_possible_".$id;$HYPCris_nbr++}
$crisprs_nbr=create_file($file_new, $file_old,$l1,$l2,$id,$CrisprBeg,$CrisprEnd,$nbsp1,$modf,1,$crisprs_nbr,@divi);
create_spFile($dir,$ResultDir,$RefSeq,$Spfile_new, $Spfile_old,0, $nbsp1);
}else{$crisprs_nbr--;}

$l1 = $line + 2;
$l2 = $j ;
$id = $rank + 1;

$Spfile_new = "Spacers_test_".$id;

# get the wrong spacer length
my $sp = getSp($file_old,$l1);

$CrisprBeg = $sp - $DRlength;
$CrisprEnd = $endPos;
my $sp_prev = $nbsp1+1;
$nbsp2 = ($l2 - $l1)/2;
if($nbsp2>=1){
if($nbsp2 >= 3){$file_new = "$inputfile"."_test_".$id;}
else{$file_new = "$inputfile"."_test_possible_".$id;$HYPCris_nbr++;}

$crisprs_nbr =create_file($file_new, $file_old,$l1,$l2,$id,$CrisprBeg,$CrisprEnd,$nbsp2,$modf,2,$crisprs_nbr,@divi);
create_spFile($dir,$ResultDir,$RefSeq,$Spfile_new, $Spfile_old, $sp_prev,$nbsp2);
}else{$crisprs_nbr--;}

$rank = $id;

# renaming the files
my $l;
if($nbsp1>=1 && $nbsp2>=1){
for(my $k=$crisprs_nbr-1; $k>=$rank; $k--){
   $l = $k + 1;
   if(-e "$dir/$RefSeq"."_Crispr_".$k){rename("$dir/$RefSeq"."_Crispr_".$k, "$dir/$RefSeq"."_Crispr_".$l);
rename("$dir/Spacers_".$k, "$dir/Spacers_".$l);
   }else{rename("$dir/$RefSeq"."_PossibleCrispr_".$k, "$dir/$RefSeq"."_PossibleCrispr_".$l);}
}
}else{

for(my $k=$crisprs_nbr-1; $k>=$rank; $k--){
   $l = $k - 1;
   if(-e "$dir/$RefSeq"."_Crispr_".$k){rename("$dir/$RefSeq"."_Crispr_".$k, "$dir/$RefSeq"."_Crispr_".$l);
rename("$dir/Spacers_".$k, "$dir/Spacers_".$l);
   }else{rename("$dir/$RefSeq"."_PossibleCrispr_".$k, "$dir/$RefSeq"."_PossibleCrispr_".$l);}
}

}

$l = $rank - 1;

if($nbsp2>=1){
if($nbsp1>=1){$l=$rank;}
unlink("$RefSeq"."_Crispr_".$rank);
if(-e "$inputfile"."_test_".$rank){rename("$inputfile"."_test_".$rank,"$RefSeq"."_Crispr_".$l);}else{rename("$inputfile"."_test_possible_".$rank,"$dir/$RefSeq"."_PossibleCrispr_".$l);}

# ___spacers files ----------
 $Spfile_new = "Spacers_test_".$rank; 
 my $Spfile = "Spacers_".$l; 
 rename("$dir/$Spfile_new", "$dir/$Spfile");
}

$l = $rank - 1;
if($nbsp1 >=1){ 
unlink("$RefSeq"."_Crispr_".$l);
if(-e "$inputfile"."_test_".$l){rename("$inputfile"."_test_".$l,"$RefSeq"."_Crispr_".$l);}
else{rename("$inputfile"."_test_possible_".$l,"$RefSeq"."_PossibleCrispr_".$l);}


# ___spacers files ----------

 $Spfile_new = "Spacers_test_".$l; 
 my  $Spfile = "Spacers_".$l; 
 rename("$dir/$Spfile_new", "$dir/$Spfile");
}

	}

  return ($rank,$crisprs_nbr,$HYPCris_nbr);
}
#--------------------
sub create_file{
 my($file_new, $file_old,$l1,$l2,$id,$CrisprBeg,$CrisprEnd,$nbspacers,$modf,$ordre,$crisprs_nbr,@divi) = @_;
 my $File_Content='';
 my @temp;
  open(FD, $file_old)  or die "Error in opening the file:$file_old!";
   while(<FD>){
      if($. <= 14){$File_Content .= $_;}

      if($. == 15){
	$File_Content .= "# Crispr Rank in the sequence: $id\n";
        $File_Content .= "# Crispr_begin_position: $CrisprBeg\t Crispr_end_position: $CrisprEnd";
      }

      if($. == 17){
	@temp = split(/:/,$_);
	$File_Content .= $temp[0];
	$File_Content .= ":";
	$File_Content .= $temp[1];
	$File_Content .= ":";
	$File_Content .= $temp[2];
	$File_Content .= ":";
        $File_Content .= " $nbspacers\n";
        $File_Content .= "#=========================================================================\n";
        $File_Content .= "Spacer_begin_position\t Spacer_length\t Spacer_sequence\n";
      }
      if ( ($.>=$l1) && ($. < $l2) ){$File_Content .= $_;}
  }
  $File_Content .= "#=========================================================================\n";


 if(($#divi == -1)||($ordre == 1)){$File_Content .="########################################\n";}

else{
    $File_Content .= "Spacers divisions:";
    foreach(@divi){$File_Content .= $_;$File_Content .=":";}
    $File_Content .="\n########################################\n";

  }

if($ordre == 2){$crisprs_nbr++;}


  close(FD);
  open WRITER,"> $file_new" or die "The file $file_new cannot be edited !\n";
  print WRITER $File_Content;
  close(WRITER);

return $crisprs_nbr;
}
#####################################
sub create_spFile{
 my ($dir,$ResultDir,$RefSeq,$Spfile_new,$Spfile_old,$sp_prev,$nbsp) = @_;
 $sp_prev = $sp_prev *2 +1;
 $nbsp = $sp_prev + $nbsp *2 -1;


 my $File_Content = "";
 open(FD, "$dir/$Spfile_old")  or die "Error in opening the file:$Spfile_old!";

   while(<FD>){
	if( ($. >= $sp_prev) && ($. <= $nbsp) ){
	  $File_Content .= $_;
	}
   }
 close(FD);
 open WRITER,"> $dir/$Spfile_new" or die "The file $Spfile_new cannot be edited !\n";
 print WRITER $File_Content;
 close(WRITER);
 
}
#------------------------------------------------------------------------------
sub getSp{
 my($file_old,$l1,) = @_;
 my (@temp, $Sp);
  
  open(FD, $file_old)  or die "Error in opening the file:$file_old!";
   while(<FD>){
      if ($. == $l1){@temp = split(/ +/,$_);}
  }
  close(FD);
  $Sp = $temp[1];
#print "taille : $Sp\n";
return $Sp;
}
#------------------------------------------------------------------------------

sub find_Sp_lines{
 my($j,$inputfile,$nb_div) = @_;
 my $i=0;
 my $line;

 open(FD, $inputfile)  or die "Error in opening the file:$inputfile!";
 while(<FD>){
     if(($_ =~ $nb_div)&&($. != $j+1)){$line= $. ;} 
 }
 close(FD);
return $line;
}
#------------------------------------------------------------------------------
sub modifySumUpfile{
 my ($dir,$ResultDir,$RefSeq,$newCrisNbr,$HYPCris_nbr) = @_;
my $inputrefseqfile = $RefSeq."_CRISPRs";
 my $tmp = "tmp";
 my $File_Content = "";
 open(FD, "$dir/$inputrefseqfile")  or die "Error in opening the file:$inputrefseqfile!";
   while(<FD>){
      if($. == 10){$File_Content .= "# Nbr_of_good_CRISPRs : $newCrisNbr\n";}
      elsif($. == 11){$File_Content .= "# Nbr_of_one_Spacer_CRISPRs : $HYPCris_nbr\n";}
      else{$File_Content .= $_;}
   }
 close(FD);
 open WRITER,"> $dir/$tmp" or die "The file $tmp cannot be edited !\n";
 print WRITER $File_Content;
 close(WRITER);
 my $rm = "rm $dir/$inputrefseqfile"; system($rm);
 rename("$dir/$tmp", "$dir/$inputrefseqfile");

}
#------------------------------------------------------------------------------
sub create_recap{
  use Date::Calc qw(:all);
  my($RefSeq, $nbrcris, $OneSpacerCris_nbr, $ResultDir) = @_;
  my $File_Content;
  my $Crispr_report = $RefSeq."_CRISPRs"; 
  my $directory = $ResultDir;
  #print $directory;
  unless(-d $directory){ mkdir $directory or die "$0: I can not create the folder $directory: $!\n" }
  my $dir = $directory."/".$RefSeq;

  if($nbrcris>0 || $OneSpacerCris_nbr>0){
      
      
  

  unless(-d $dir){ mkdir $dir or die "$0: I can not create the folder $dir: $!\n" }
  
  my($year,$month,$day, $hour,$min,$sec) = Today_and_Now(1);
    $File_Content .= "#=======================================\n";
  $File_Content .= "# \n";
  $File_Content .= "# Sequence: $RefSeq \n";
  $File_Content .= "# Number of CRISPR-like structures : $nbrcris\n";
  $File_Content .= "# Number of questionable structures : $OneSpacerCris_nbr\n";
  $File_Content .= "# \n";
  $File_Content .= "#=======================================\n";
  #chdir($dir);
  open WRITER,"> $dir/$Crispr_report" or die "The file $Crispr_report cannot be edited !\n";
  print WRITER $File_Content;
  close WRITER;
  }#else{
    #  my $dirt="tmp";
   #unless(-d $dirt){ mkdir $dirt or die "$0: I can not create the folder $dirt: $!\n" }
    #  chdir($dirt);
  #}
}
#------------------------------------------------------------------------------
# copy spacers file of crispr_$refseq_$i
# this file will be edited by user (button)
sub copy_spacers{

  my($SpacersFile,$dir,$id) = @_; 
  $dir = $dir."/Spacers_$id";
#print "j'ai copiÈ $SpacersFile dans $dir";
  my $cop = "cp -R $SpacersFile $dir";
  system($cop);
}
#------------------------------------------------------------------------------
# fill in the file crispr_$refseq_$i
# this file will be used for the database storage
sub fill_in_crisprfile
{
  use Date::Calc qw(:all);

  my(	$simDRs,$RefSeq,$ResultDir,$id,$CrisprBeg,$CrisprEnd,
  	$DRseq,$nbspacers,$SpacersFile,$RefFalsSpacers,%spacersPos) = @_; 


  my ($Crispr_file,$File_Content,$Bpos,$i);
  my $DRlength = length($DRseq);
  my($year,$month,$day, $hour,$min,$sec) = Today_and_Now(1);
  my $directory = $ResultDir;
  if(-e "$directory"){}else{mkdir($directory);}
  my $dir = $directory."/".$RefSeq;
  if(-e "$dir"){}else{mkdir($dir);}

  if( ($nbspacers >= 2) && ($simDRs ==1) )
  {
  	$Crispr_file = "$dir/$RefSeq"."_Crispr_".$id;
  }
  else
  {
  	$Crispr_file = "$dir/$RefSeq"."_PossibleCrispr_".$id;
  }
  $nbspacers ++;
  $File_Content  = "########################################\n";
  $File_Content .= "# Program: CRISPR Finder\n";
  $File_Content .= "# Author: Ibtissem GRISSA\n";
  $File_Content .= "# Rundate (GMT): $day/$month/$year $hour:$min:$sec\n";
  $File_Content .= "# Report_file: $Crispr_file\n";
  $File_Content .= "########################################\n";
  $File_Content .= "#=======================================\n";
  $File_Content .= "# \n";
  $File_Content .= "# Sequence: $RefSeq \n";
  $File_Content .= "# Description: ".$seq->desc()."\n";
  $File_Content .= "# Length: ".$seq->length()."\n";
  $File_Content .= "# Id: ".$seq->id()."\n";
  $File_Content .= "#\n";
  $File_Content .= "#=========================================================================\n";
  $File_Content .= "# Crispr Rank in the sequence: $id\n";
  $File_Content .= "# Crispr_begin_position: $CrisprBeg\t Crispr_end_position: $CrisprEnd\n";
  $File_Content .= "# DR: $DRseq\t DR_length: $DRlength\t Number_of_spacers: $nbspacers\n";
  $File_Content .= "#=========================================================================\n";
  $File_Content .= "Spacer_begin_position\t Spacer_length\t Spacer_sequence\n";

  copy_spacers($SpacersFile, $dir,$id);
  open(FD, $SpacersFile)  or die "Error in opening the input file : $SpacersFile!";
  my @lines=<FD>;
  close(FD);
  $i=1;
  foreach my $k (reverse (sort { $b <=>$a } keys %spacersPos) )
  {
	$Bpos = " "x(20 - length($k));
	$File_Content .= "$Bpos $k";
	$Bpos = ' 'x(12 - length($spacersPos{$k}));
	$File_Content .= "\t $Bpos $spacersPos{$k}";
	$File_Content .= "\t $lines[$i]\n";
	$i +=2;
  }
  $File_Content .= "#=========================================================================\n";
  my @FalsSpacers = @$RefFalsSpacers;
  if($#FalsSpacers > -1)
  {
    $File_Content .= "Spacers divisions:";
    foreach(@FalsSpacers){$File_Content .= $_;$File_Content .=":";}
    $File_Content .="\n";
  }
  $File_Content .="########################################\n";
  open WRITER,"> $Crispr_file" or die "The file $Crispr_file cannot be edited !\n";
  print WRITER $File_Content;
  close WRITER;
  return $Crispr_file;
}

#------------------------------------------------------------------------------
# compute an error score for a DR matches in the crispr file
# used in the comparison of two DRs
sub compute_Crispr_score
{
  my($crisprfile,$DRlength) = @_;
  my $score = 0;
  my(@lines, @temp, $penal, $TotERR,$nb);
  my $troncated = 0;  
  $TotERR =0; $nb=0;
#print "<br>";
  open(FD, $crisprfile)  or die "Error in opening the input file : $crisprfile!";
  @lines=<FD>;
  close(FD);
  my $count = 19;
  if($#lines <= 5){$score = 100000;}	#added in 15/11/06	
  for($count=12; $count < scalar(@lines)-1; $count++)
  {
     $nb++;
     next if $lines[$count] =~ /^#/;
     next if $lines[$count] =~ /^\s+$/;
     $lines[$count] = trim($lines[$count]);
     next if $lines[$count] =~ /^Start/;
     @temp = split(/ +/, $lines[$count]);
#     if($temp[0]=~""){shift @temp;}
     
     if($temp[$#temp-1] ne "."){
	$TotERR = $TotERR + $temp[$#temp-1];
	$penal = 1 + $temp[$#temp-1]/$DRlength;
#print "$count penal: $penal <br>";
	$score = $score + $penal;
	if($count == 19){$troncated = $penal;}
	if($count == scalar(@lines)-4){
		if($penal > $troncated){
			$troncated = $penal;
		}
	}
     }
  }
  #print "tronc : $troncated <br>";
  $score = $score - $troncated;
  $TotERR = $TotERR/($nb*$DRlength);
  #print "<br>$crisprfile score : $score __ toterr : $TotERR <br>";
  #$ans = <STDIN>;
  return ($score,$TotERR);
}
#------------------------------------------------------------------------------
# study a given seq_v according to a specific candidate DR
# 26/10/06 use of clustalw to align spacers
sub Find_theCrispr
{
   my($indexname,$DR,$count,$seqbeg,$seqend,$crisprfile)= @_;
   my(@spacers, $goodcrispr,$CrisBeg,$CrisEnd,$refFalsSpacers,$simDRs);
   my $DRlength = length($DR);

   ($simDRs,$refFalsSpacers, @spacers) = definespacers($DRlength, $crisprfile,$indexname);
   my @FalsSpacers = @$refFalsSpacers;
   if($#FalsSpacers >=0){for(my $h =0; $h<= $#FalsSpacers; $h++){$FalsSpacers[$h] += $seqbeg;}}
   my %spacersH = trans_struct_hach($seqbeg,@spacers);


   $goodcrispr = 0;
   my $o = $#spacers;


   if($#spacers >= 0)
   {
	if($#spacers == 0)
	{
		$goodcrispr = check_short_crispr($DR, "spacers".$indexname);
	}
	else
	{
		$goodcrispr = checkspacersAlign("spacers".$indexname);
	}
	if($goodcrispr)
	{
		$CrisBeg=0; $CrisEnd=0;
		$CrisBeg = $seqbeg + $spacers[0]->Pos1 - $DRlength;
		$CrisEnd = $seqbeg + $spacers[$#spacers]->Pos1 + $spacers[$#spacers]->Length +$DRlength -1;
	}
  }
  my $r=$#FalsSpacers;
  
  return($goodcrispr,$simDRs,$CrisBeg,$CrisEnd,\%spacersH,$#spacers,\@FalsSpacers);
}

#------------------------------------------------------------------------------
# transform the old @spacers structure to a hachage containing as keys the
# begin position of a spacer (according to the whole sequence) and as values
# the corresponding length
sub trans_struct_hach
{
   my($seqbeg,@structspacers) = @_;
   my %spacers = ();
   my ($i,$pos, $len);
   for($i=0; $i<= $#structspacers; $i++)
   {
	$pos = $structspacers[$i]-> Pos1 + $seqbeg;
	$len = $structspacers[$i]-> Length;
	$spacers{$pos} = $len;
   }
   return %spacers;
}

#------------------------------------------------------------------------------
# this function is used to check, in the case of a unique spacer, the identity 
# between the DR and the CRISPR 
# use of program needle of emboss to do the alignment

sub check_short_crispr
{
  my($DR, $spacerfile) = @_;
  my($needleoptions, $goodcrispr);

  # edit DR file
  open WRITER,"> DR" or die "The file cannot be edited !\n";
  print WRITER "> DR\n";
  print WRITER $DR;
  print WRITER "\n";
  close WRITER;

  # needle program
  $needleoptions = "needle -auto -asequence DR -bsequence $spacerfile -gapopen 10.0 -gapextend 0.5 -outfile alignDR_Spacer.needle";
  makesystemcall($needleoptions);
  $goodcrispr = similarity_needle();
#print "goodneedle = $goodcrispr <br>";
  return $goodcrispr;
}
#------------------------------------------------------------------------------
# analyze needle file to conclude if repetitivity is significant
sub similarity_needle{

  open(FD, "alignDR_Spacer.needle")  or die "Error in opening the input file : alignDR_Spacer.needle!";
  # find the gap value
  while(<FD>){
	if($_ =~ "Gaps:" ){
	  my @temp = split(/ +/, $_);
	  my $score = $temp[3];
	  $score =~ s/\(//;
	  $score =~ s/%\)//;
	  if($score > 50.0){return 1;}else{return 0;}
	}
  }
  close(FD);
}

#------------------------------------------------------------------------------
# checks if the DR is already a candidate 
sub check_DR{
  my ($DR, @DR_cand) = @_;
  my $i=0;
  my $stop = 0;
  
  while(($i <= $#DR_cand) && ($stop==0) ){
	if($DR eq $DR_cand[$i]){$stop = 1;};
	$i++;
  }
return $stop;
}
#------------------------------------------------------------------------------
# counts the number of occurences of the input DR in the analyzed genome 
sub DR_occ{
  my ($DR, @rep) = @_;
  my $DRocc = 0;
  my $count=0;
  for($count=0; $count<=$#rep; $count++){
	if($DR eq $rep[$count]->DRseq){
		$DRocc++;
	}
  }
  return $DRocc;
}
#------------------------------------------------------------------------------
# counts the number of occurences of the input DR in the analyzed genome (the reverse complement is counted also)
sub DR_occ_rev
{
  my ($DR, @rep) = @_;
  my $DRocc = 0;
  my $count=0;
  my($seqDR, $revDR);

  $seqDR = Bio::Seq->new(-seq => $DR);
  $revDR = $seqDR->revcom();
  for($count=0; $count<=$#rep; $count++){
	if(($DR eq $rep[$count]->DRseq) || ($revDR->seq()eq$rep[$count]->DRseq) ){
		$DRocc++;
	}
  }
  return $DRocc;
}
#######################################################################
# create the file DRconsenus_i
sub print_DR
{
  my($count, $seqbeg, $seqend, @DR_array) = @_;
  if($#DR_array <1 ){
     my $file_DRconsensus = "DRconsensus_" . $count;
     open WRITER,"> $file_DRconsensus" or die "The file $file_DRconsensus cannot be edited !\n";
     print WRITER "> Hypothetical DR for the seq [$seqbeg, $seqend]\n";
     print WRITER $DR_array[0];
     print WRITER "\n";
     close(WRITER);
     return 0;
  }else
  {
     my $i=0;
     while($i <= $#DR_array){	
	my $file_DRconsensus = "DRconsensus_" . $count."_$i";
	open WRITER,"> $file_DRconsensus" or die "The file $file_DRconsensus cannot be edited !\n";
        print WRITER "> Hypothetical DR for the seq [$seqbeg, $seqend]\n";
        print WRITER $DR_array[$i];
        print WRITER "\n";
        close(WRITER);
        $i++;
     }
     $i = $i-1;
     return $i;
  }
}
#------------------------------------------------------------------------------
#find clusters and store their start and end positions
# @tabseq = find_clusters(@rep);
sub find_clusters		# find all possible clusters (even the shortest)
{
 my @repetitions = @_;
 my @tabseq = ();
 my $count = 0;
 my $nb_clust = 0;	# nb_clust = 2*nbr_clusters

#print "nbr de rep = $#repetitions --\n";
$tabseq[$nb_clust] = $repetitions[$count]->Pos1;	# start position of the first cluster
  while($count < $#repetitions)
  {
    if(!compare_clusters($repetitions[$count]->Pos1,$repetitions[$count+1]->Pos1))
    {
     # a new cluster is found
     if($repetitions[$count]->Pos1 !=$tabseq[$nb_clust])
     {
	$nb_clust++; 
	$tabseq[$nb_clust] = $repetitions[$count]->Pos1;
      }
      else
      {
	$nb_clust++; 
	$tabseq[$nb_clust] = $repetitions[$count]->Pos1 + 60;

      }

     $nb_clust++; 
     $tabseq[$nb_clust] = $repetitions[$count+1]->Pos1;
	
     }
     $count++;
   }


   if($#tabseq % 2 == 0)
   {
		if( $repetitions[$count]->Pos1 !=$tabseq[$nb_clust] )
		{
			$nb_clust++; 
			$tabseq[$nb_clust] = $repetitions[$count]->Pos1;
		}
		else
		{
			$nb_clust++; 
			$tabseq[$nb_clust] = $repetitions[$count]->Pos1 + 60;
		}
   }
   return @tabseq;
}

#------------------------------------------------------------------------------
sub compare_clusters
{
  my($element1,$element2) = @_;
  if( ($element2 >= $element1 - 1500) && ($element2 <= $element1+1500) )
  {
    return 1;
  }else
  {
    return 0;
  }
}

#------------------------------------------------------------------------------
sub checkdbfile
{
  my($inputfile,$prjfile) = @_;
  my $dbfile = '';

  if(-e $prjfile)
  {
    unless(open(PRJFILEPTR,$prjfile))
    {
      return 0;
    }
    while(my $line = <PRJFILEPTR>)
    {
      if($line =~ /^dbfile=(\S+) (\d+)/)
      {
        if($dbfile eq '')
        {
          $dbfile = $1;
          my $dbfilesize = $2;
          if($dbfile eq $inputfile)
          {
            if(my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,
                   $atime,$mtime,$ctime,$blksize,$blocks) = stat($dbfile))
            {
              if($size == $dbfilesize)
              {
                close PRJFILEPTR;
                return 1;
              }
            }
          }
        }
      }
      close PRJFILEPTR;
      return 0;
    }
  }
  return 0;
}
#------------------------------------------------------------------------------
sub makesystemcall
{
  my ($argstring) = @_;
  my @args = split(' ',$argstring);
  my($retcode) = system($argstring);
  $retcode = $? >> 8;
  if($retcode ne 0)
  {
    print STDERR "failure: \"$argstring\", errorcode $?\n";
    exit 1;
  }
  #print STDERR "# $argstring\n"; #skip verbose mode
}

sub callmkvtree
{
  my ($inputfile,$indexname) = @_;
  if(not (checkdbfile($inputfile,"${indexname}.prj")))
  {
    #makesystemcall("./mkvtree -db $inputfile " .
    makesystemcall("mkvtree -db $inputfile " .   # LK
                   "-dna -pl -lcp -suf -tis -ois -bwt -bck -sti1");
  }
}
#------------------------------------------------------------------------------
sub analyzeCFargs
{
  my ($program,$argvref,$argcount) = @_;
  my @vmatchoptions = ();
  my %notsupported = 
  (
    "-r" => 1,
    "-c" => 1,
    "-hrate" => 1,
    "-erate" => 1,
    "-o" => 1,
    "-b" => 1,
    "-warn" => 1,
    "-iw" => 1,
    "-mem" => 1
  );

  my @argv = @$argvref;
  my $argnum;
  my $stringoption = 0;
  my $linewidth = 0;
  my $doiub = 0;
  my $bestoption = 0;
  my $allmaxoption = 0;
  for($argnum=0; $argnum<$argcount-1; $argnum++)
  {
    my $currentarg = $argv[$argnum];
    if($currentarg eq '-f')
    {
      push(@vmatchoptions,"-d");
    } elsif($currentarg eq '-p')
    {
      push(@vmatchoptions,"-p");
    } elsif($currentarg eq '-l' || 
            $currentarg eq '-seedsize' ||
            $currentarg eq '-best')
    {
      if($currentarg eq '-seedsize')
      {
        push(@vmatchoptions,"-seedlength");
      } else
      {
        push(@vmatchoptions,$currentarg);
      }
      $argnum++;
      if($argnum >= $argcount-1 || $argv[$argnum] =~ m/^\-/)
      {
        print STDERR "$program: missing argument for option \"$currentarg\"\n";
        exit 1;
      }
      if($currentarg eq '-best')
      {
        $bestoption = 1;
      }
      push(@vmatchoptions,$argv[$argnum]);
    } elsif($currentarg eq '-lw')
    {
      $argnum++;
      if($argnum >= $argcount-1 || $argv[$argnum] =~ m/^\-/)
      {
        print STDERR "$program: missing argument for option \"$currentarg\"\n";
        exit 1;
      }
      $linewidth = $argv[$argnum];
      if($linewidth <= 0)
      {
        print STDERR "$program: illegal argument \"$linewidth\" ";
        print STDERR "to option \"-lw\"\n";
        exit 1;
      }
    } elsif($currentarg eq '-h' || $currentarg eq '-e')
    {
      push(@vmatchoptions,$currentarg);
      $argnum++;
      if($argnum >= $argcount-1 || $argv[$argnum] =~ m/^-/)
      {
        push(@vmatchoptions,"4");
      } else
      {
        push(@vmatchoptions,$argv[$argnum]);
      }
    } elsif($currentarg eq '-allmax')
    {
      $allmaxoption = 1;
      push(@vmatchoptions,"-allmax");
    } elsif($currentarg eq '-s')
    {
      $stringoption = 1;
    } elsif($currentarg eq '-iub')
    {
      $doiub = 1;
    } elsif($currentarg eq '-nodistance')
    {
      push(@vmatchoptions,"-nodist");
    } elsif($currentarg eq '-noevalue' || $currentarg eq '-i')
    {
      push(@vmatchoptions,$currentarg);
    } else
    {
      if(exists $notsupported{$currentarg})
      {
        print STDERR "$program: CF option \"$currentarg\" is ";
        print STDERR "not supported\n";
      } else
      {
        print STDERR "$program: illegal option \"$currentarg\"\n";
      }
      exit 1;
    }
  }
  if($argnum == $argcount-1)
  {
    if($argv[$argnum] =~ m/^-/)
    {
      print STDERR "$program: last argument must be filename, ";
      print STDERR "not beginning with \"-\"\n";
      exit 1;
    } 
  } elsif($argnum > $argcount-1)
  {
    print STDERR "$program: missing last argument\n";
    exit 1;
  } elsif($argnum < $argcount-1)
  {
    print STDERR "$program: superfluous argument $argv[$argnum]\n";
    exit 1;
  }
  if(not @vmatchoptions)
  {
    print STDERR "$program: at least one option is required\n";
    exit 1;
  }
  if((not $bestoption) && (not $allmaxoption))
  {
    push(@vmatchoptions,"-best");
    push(@vmatchoptions,"50");
  }
  if($stringoption)
  {
    push(@vmatchoptions,"-s");
    if($linewidth > 0)
    {
      push(@vmatchoptions,$linewidth);
    }
    if($doiub)
    {
      push(@vmatchoptions,"abbreviub");
    }
  }
  push(@vmatchoptions,"-noscore");
  push(@vmatchoptions,"-noidentity");
  push(@vmatchoptions,"-absolute");
  return \@vmatchoptions;
}

#------------------------------------------------------------------------------

#find a cluster by starting search from position index
#retuns the index end of the cluster and creates  else return 0
sub find_cluster	# version 25/04/2006
{
 my @arg = @_;
 my(@repetitions,@clusters_def);        #@clusters_def[0] is the index from which we begin clusters search, 
 my $i;                                    #the other values of the array (when they exist) are resp. first_pos, length of hypothetical crisprs loci
 for( $i = 1; $i<=  $arg[0]+1; $i++){$repetitions[$i-1] = $arg[$i];}
 my $j=0;
 for($i = $arg[0]+2; $i<= $#arg; $i++)
 {
	$clusters_def[$j] = $arg[$i];
	$j++;
 }

 my($pos_min, $pos_end, @DR_min, $nbr_clusters, $file_DRconsensus, $ref);
 my $count = $clusters_def[0];
 my $test = 0;
 my $n_DR = 0;
 my $k=0;

 while( ($count < $#repetitions) && ($test==0) )
 {

	if(!compare_clusters($repetitions[$count]->Pos1,$repetitions[$count+1]->Pos1) )
	{
	   $count++;

	}
	else
	{
	    $test = 1;
	    $ref = $count;

	}
  }
  if($test == 1)   # a new cluster is found
  {
		$pos_min = $repetitions[$ref]->Pos1;
		$DR_min[0] = $count;
	
		while( (compare_clusters($repetitions[$ref]->Pos1,$repetitions[$count+1]->Pos1)) )
		{
	    	if($repetitions[$DR_min[0]]->Length > $repetitions[$count+1]->Length)
	    	{
				@DR_min = ();
				$DR_min[0] = $count+1;
	    	}
	    	# case we have different DRs of the same size
	    	if($repetitions[$DR_min[0]]->Length == $repetitions[$count+1]->Length)
	    	{
				if(!($repetitions[$DR_min[0]]->DRseq eq $repetitions[$count+1]->DRseq))
				{ 
		   		while( $k <= $#DR_min )
		   		{
						if(!($repetitions[$DR_min[$k]]->DRseq eq $repetitions[$count+1]->DRseq)){push(@DR_min, $count+1);}   #add the new DR to DRs array
						$k++;
		   		}
				}
	    	}
	    	$count++;
	    	if($count == $#repetitions){last;}
		}
	
		# fill in the table of hypothetical crisprs loci
		$pos_end =  $repetitions[$count]->Pos1 + 2*$repetitions[$DR_min[0]]->Length+1;
		push(@clusters_def, $pos_min);  # add the first position
		push(@clusters_def, $pos_end);      # add the minimum sequence length
	
		# create DR consensus file
		$nbr_clusters = $#clusters_def/2 ;
		$file_DRconsensus = "DRconsensus_" . $nbr_clusters;
		open WRITER,"> $file_DRconsensus" or die "The file $file_DRconsensus cannot be edited !\n";

		# case we have many possible DRs of the same length
		if($#DR_min >0)
		{
			print"\n ======There are $#DR_min+1 possible DRs of the same length in the sequence $file_DRconsensus =====\n\n";

			for($k=0; $k<=$#DR_min; $k++)	
			{
				print WRITER "> Hypothetical Direct repeat for the region around $pos_min\n";
				print WRITER $repetitions[$DR_min[$k]]->DRseq;
				print WRITER "\n";
			}
			close(WRITER);
			
			callmkvtree($file_DRconsensus,$file_DRconsensus);
			my $taille = length($repetitions[$DR_min[0]]->DRseq) -2;  #case we have two differences
			my $a="vmatch -l ".$taille . " -noevalue -noscore -supermax -noidentity -s leftseq ".$file_DRconsensus." > newtemp";
			print $a;
			# makesystemcall("./vmatch -l ".$taille . " -noevalue -noscore -supermax -noidentity -s leftseq ".$file_DRconsensus." > newtemp"); 
			print("vmatch -l ".$taille . " -noevalue -noscore -supermax -noidentity -s leftseq ".$file_DRconsensus." > newtemp");
			makesystemcall("vmatch -l ".$taille . " -noevalue -noscore -supermax -noidentity -s leftseq ".$file_DRconsensus." > newtemp");  #LK
			  
			my @temp;
			open(FD, "newtemp")  or die "Error in opening the input file : $!";
 			while( <FD> ) 
			{
            if($. == 3)
		   	{
					open WRITER,"> $file_DRconsensus" or die "The file $file_DRconsensus cannot be edited !\n";
					print WRITER "> Hypothetical Direct repeat for the region around $pos_min\n";
					print WRITER $_;
					print WRITER "\n";
					close(WRITER);
		   	}
  			}
			close(FD);
			unlink("newtemp") or die " can't delete the temporary file newtemp : !\n";
			callmkvtree($file_DRconsensus,$file_DRconsensus);
		}
		else
		{
			print WRITER "> Hypothetical Direct repeat for the region around $pos_min\n";
			print WRITER $repetitions[$DR_min[0]]->DRseq;
			print WRITER "\n";
			close(WRITER);
		}
			#put the next index for the following cluster search
			$clusters_def[0] = $count; 
  }
   else{$clusters_def[0] = -1;}
   return @clusters_def;
}
#-------------------------------------------------------------------------------
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
#------------------------------------------------------------------------------
sub printhelpall
{
print <<HEREDOC;
Command line: CRISPRFinder-v3.pl <filename>
Your result folder will be in the same directory

Command line: CRISPRFinder-v3.pl <filename> <directory>
Your result folder will be in the chosen directory

Sequences should be in fasta format,
The program requires mkvtree and vmatch installed.
vmatch: copyright Stefan Kurtz 2000-2005, kurtz\@zbh.uni-hamburg.de
Dependencies:
- perl,
- perl's modules: Class::Struct, Bio::Tools::Run::Alignment::Clustalw , Date::Calc, File::Copy, Bio::Seq and Bio::SeqIO, 
- EMBOSS 5.0.0 or upper
HEREDOC
}

sub printversion
{
  my $programname = shift @_;
  
  print "\nThis is $programname, version 3,\n";
  print " a perl script to identify CRISPRs in DNA sequences\n";
  
}

sub printhelpbasic
{
  my $programname = shift @_;
  print STDERR "Usage: $programname [options] <filename> <Results_directory>\n";
  print "type -v to see the current version\n";
  print "type -help for help\n";
}


