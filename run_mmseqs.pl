#!/usr/bin/perl

# Run MMSeqs search for stitched samples against the Refseq complete Protein database

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

$USAGE = q/USAGE:
perl run_mmseqs.pl    -db <\/path\/to\/MMseqsReferenceDB> 
                      Complete path to the MMseqs DB that will 
                      be used for the search (E.g an MMseqs DB made from
                      NCBI Refseqs Complete OR Uniref90 etc.)
                     
                      -p, --parallel <Int number>
                      Number of parallel processes to run

                      -o, --out <\/path\/to\/store\/results>    
                      Name of the output folder where the results and
                      intermediate files will be stored 
                      Default: \"\.\"
                   
                      -m, --meganpath <\/path\/to\/Megan\/Installation>
                      Path to the Megan installation folder to run Megan.
                      This folder should contain the \/tools directory with 
                      the megan scripts blast2rma and rma2info
                      This folder should also contain the mapping files provided by 
                      Megan
                      If not specified, only the mmseqs steps of the pipeline 
                      will be run
                     
  
/;


use Getopt::Long;
use Data::Dumper;
use Parallel::ForkManager;
use Sys::CPU;


GetOptions (
'o|out=s' => \$outpath,
'db=s' => \$db,
'p|parallel:i' => \$parallel,
'm|meganpath=s' => \$meganpath
) or die $USAGE;

die $USAGE if !$db;

@files = @ARGV;

print scalar @files,"\n";

    if (!$outpath)
    {
    $outpath=".";
    }
    else
    {
    system ("mkdir -p $outpath");
    }

my $cpu_count=1;
#if the option is set
    if(defined($parallel))
    {
        #option is set but with no value then use the max number of proccessors
        if($parallel == 0)
        {
            #load this module dynamically
            eval("use Sys::CPU;");
            $cpu_count=Sys::CPU::cpu_count();
        }
        else
        {
            $cpu_count=$parallel;
        }
    }

$base_cmd="mmseqs";
$dbIndex_cmd = $base_cmd." createindex $db tmp";

print "Total Cpus usd:$cpu_count\nRunning mmseqs index on $db ... $dbIndex_cmd\n";

#system("$dbIndex_cmd > /dev/null");
#system("mmseqs createindex $db tmp");
#system("mmseqs createindex /home/dhwani/databases/mmseqsRefSeqCompleteDB tmp");

$pm = new Parallel::ForkManager($cpu_count);

$pm->run_on_finish( sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
    
    my $q = $data_structure_reference->{input};
    
    print "Finished running process $pid for $q, $exit_code, $exit_signal\n";
#     $results{$q} = $data_structure_reference->{result};
});

    foreach my $file (@files) 
    {
    ($path,$fname) = $file =~ m|^(.*[/\\])([^/\\]+?)$|;

        if ($fname =~ /_/)
        {
        @temp=split("_",$fname);
        $sample_tag = $temp[0];
        }
        else
        {
        $sample_tag = $fname;
        $sample_tag =~ s/\.fastq//;
        }

    $querydb=$base_cmd."-".$sample_tag."queryDB";
    $result_db=$base_cmd."-".$sample_tag."resultDB";
    $m8_outfile=$base_cmd."-".$sample_tag."-s1.m8";
    $rma_outfile=$base_cmd."-".$sample_tag."-top25-tax-interpro.rma";
    $megantax_out=$base_cmd."-".$sample_tag."-Megan_Taxonomy-assignments.txt";
    $meganfunc_out=$base_cmd."-".$sample_tag."-Megan_Interpro-assignments.txt";
    $jobfile=$base_cmd."-".$sample_tag."-jobfile.sh";
    $logfile=$base_cmd."-".$sample_tag."-logfile.txt";

    open (OUT, ">$jobfile");

    print OUT "mmseqs createdb $file $outpath/$querydb\n";
    print OUT "mmseqs search $outpath/$querydb $db $outpath/$result_db tmp --db-load-mode 3 --threads 2 --max-seqs 25 -s 1 -a -e 1e-5\n";
    print OUT "mmseqs convertalis $outpath/$querydb $db $outpath/$result_db $outpath/$m8_outfile --db-load-mode 3\n";
        if ($meganpath)
        {
        print OUT "\#MEGAN will be run!\n";
        print OUT "$meganpath/tools/blast2rma -i $outpath/$m8_outfile -f BlastTab -a2t $meganpath/prot_acc2tax-Nov2018X1.abin -a2interpro2go $meganpath/acc2interpro-June2018X.abin -o $outpath/$rma_outfile\n";
        print OUT "$meganpath/tools/rma2info -i $outpath/$rma_outfile --names -r G -v -r2c Taxonomy > $outpath/$megantax_out\n";
        print OUT "$meganpath/tools/rma2info -i $outpath/$rma_outfile --names -r G -v -r2c INTERPRO2GO > $outpath/$meganfunc_out\n";
        }
        else
        {
        print "MEGAN path not provided! Not running MEGAN\n";
        } 
        
    print "\n";

    system("chmod 755 $jobfile");

    my $pid = $pm->start and next; # do the fork

    system("sh $jobfile > /dev/null 2> $logfile");
    
    my @pids = $pm->running_procs;

#     print "\n\n ----> $pid <---- @pids\n\n";

    # system("mmseqs createdb $file $outpath/$querydb");
    # system("mmseqs search $outpath/$querydb $db $outpath/$result_db tmp --db-load-mode 2 --threads 1 --max-seqs 25 -s 1 -a -e 1e-5");
    # system("mmseqs convertalis $outpath/$querydb $db $outpath/$result_db $outpath/$m8_outfile");
    # 
    #     if($meganpath)
    # 	{
    # 	system ("$meganpath/tools/blast2rma -i $outpath/$m8_outfile -f BlastTab -a2t $meganpath/prot_acc2tax-Nov2018X1.abin -a2interpro2go $meganpath/acc2interpro-June2018X.abin -o $outpath/$rma_outfile");
    # 	system ("$meganpath/tools/rma2info -i $outpath/$rma_outfile --names -r G -v -r2c Taxonomy > $outpath/$megantax_out");
    # 	system ("$meganpath/tools/rma2info -i $outpath/$rma_outfile --names -r G -v -r2c INTERPRO2GO > $outpath/$meganfunc_out");
    # 	}
    
    $pm->finish(0, { result => \@pids, input => $file });; # do the exit in the child process
    }

$pm->wait_all_children;
# print Dumper \%results;
