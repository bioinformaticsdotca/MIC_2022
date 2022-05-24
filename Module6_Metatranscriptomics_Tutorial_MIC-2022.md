  
    
<h1 dir="auto"><a id="user-content-metapro-metatranscriptomics-practical-lab" class="anchor" aria-hidden="true" href="#metapro-metatranscriptomics-practical-lab"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-metapro-metatranscriptomics-practical-lab" href="#metapro-metatranscriptomics-practical-lab"></a>MetaPro Metatranscriptomics Practical Lab</h1>
<p dir="auto"><strong>This work is licensed under a <a href="https://creativecommons.org/licenses/by-sa/4.0/" rel="nofollow">Creative Commons Attribution-ShareAlike 4.0 International</a>. This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.</strong></p>
<p dir="auto">This tutorial was developed by Billy Taj (<a href="mailto:billy.taj@sickkids.ca">billy.taj@sickkids.ca</a>), Mobolaji Adeolu (<a href="mailto:adeolum@mcmaster.ca">adeolum@mcmaster.ca</a>), John Parkinson (<a href="mailto:john.parkinson@utoronto.ca">john.parkinson@utoronto.ca</a>) &amp; Xuejian Xiong (<a href="mailto:xuejian@sickkids.ca">xuejian@sickkids.ca</a>), and updated by Ana Popovic for the IMPACTT/IMC Bioinformatics Workshop, 2022.</p>
<h2 dir="auto"><a id="user-content-overview" class="anchor" aria-hidden="true" href="#overview"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-overview" href="#overview"></a>Overview</h2>
<p dir="auto">This tutorial will take you through the MetaPro pipeline for processing and annotation of metatranscriptomic data (<a href="https://doi.org/10.1101/2021.02.23.432558" rel="nofollow">https://doi.org/10.1101/2021.02.23.432558</a>). The pipeline was developed by the Parkinson lab, and consists of the following steps:</p>
<ol dir="auto">
<li>Remove adapter sequences (added during library preparation for sequencing) and low-quality reads.</li>
<li>Remove duplicate reads to reduce processing time for downstream steps.</li>
<li>Remove vector contamination (reads derived from cloning vectors, spike-ins, and primers).</li>
<li>Remove host reads (if exploring a host-derived microbiome).</li>
<li>Remove rRNA sequences, which typically dominate metatranscriptomic datasets.</li>
<li>Repopulate duplicate reads (removed in step 2).</li>
<li>Assemble reads into contigs, and predict ORFs.</li>
<li>Annotate ORFs to known genes and proteins.</li>
<li>Classify reads to known taxonomic groups, and visualize the taxonomic composition of your dataset.</li>
<li>Predict enzyme function.</li>
<li>Generate output files, including gene annotations with normalized gene expression values, and enzyme predictions.</li>
<li>Visualize metabolic pathway activity, as mapped onto KEGG-defined pathways, within Cytoscape.</li>
</ol>
<p dir="auto">The MetaPro metatranscriptomics pipeline includes a series of Python scripts that handle the orchestrated invocation of existing bioinformatics tools, read-conflict resolution, file format conversion, and output parsing. We will go through these steps to illustrate the complexity of the process and the underlying tools and scripts. The pipeline is designed to run all steps automatically in succession.  However, a tutorial mode is available to allow running of single steps, which we will employ in this tutorial for the purpose of learning.</p>
<p dir="auto">New, faster, and/or more accurate tools are being developed all the time, and it is worth bearing in mind that any pipeline needs to be flexible to incorporate these tools as they are adopted as standards by the community. For example, over the past two years, our lab has transitioned from using Trimmomatic to AdapterRemoval, and from BLAST to DIAMOND.</p>
<p dir="auto">To illustrate the process, we are going to use sequenced reads generated from the colon contents of a mouse. The original dataset contains 25 million 150 bp single-end reads. Rather than use the entire set, which might take several days to process on a desktop, the tutorial will take you through processing a subset of 100,000 single-ended reads. Paired-end reads can also be used, and are often preferred because they can improve annotation quality when there is enough overlap in read pairs to improve the effective average read length. Working with paired-end data involves an additional data processing step (merging of overlapping reads) and produces more files (files for merged/singleton reads, forward reads and reverse reads), but the structure of the remaining steps is the same as that described below.</p>
<p dir="auto"><strong>Note:</strong>
The purpose of this tutorial is to demonstrate MetaPro's various steps. The pipeline is fully capable of running all steps without the intervention of the user, beyond an initial call to the program. If you wish to use automated MetaPro pipeline, do not use the --tutorial option.</p>
<h2 dir="auto"><a id="user-content-preliminaries" class="anchor" aria-hidden="true" href="#preliminaries"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-preliminaries" href="#preliminaries"></a>Preliminaries</h2>
<h3 dir="auto"><a id="user-content-launch-the-metapro-pipeline" class="anchor" aria-hidden="true" href="#launch-the-metapro-pipeline"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-launch-the-metapro-pipeline" href="#launch-the-metapro-pipeline"></a>i. Launch the MetaPro pipeline</h3>
<p dir="auto">MetaPro operates in a containerized environment running Ubuntu Linux 18.04. The pipeline and its dependent programs exist within a Docker container image, and may be accessed using the Docker or Singularity tools. For the purposes of this tutorial, Docker and Singlarity have both already been installed, and the MetaPro pipeline has been downloaded to your environment.</p>
<p dir="auto">Docker and Singularity maintain different access modes to use their containers.</p>
<ol dir="auto">
<li>Scripted-mode: where the user can call software within the container, to run it.</li>
<li>Interactive-mode: where the user can enter into the container and use it like an operating system.</li>
</ol>
<p dir="auto">In this tutorial, we will be using Singularity in interactive mode to run MetaPro commands. The command to download and create the MetaPro image <em>(below)</em> has already been run for you:</p>
<p dir="auto">    singularity pull docker://parkinsonlab/metapro:develop   <strong>[DO NOT RUN]</strong></p>
<p dir="auto"><br><br>Navigate to your workspace and create a new tutorial folder (we will be using the <em>absolute</em> path):</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="cd /media/data/workspace  
mkdir metapro_tutorial/  
cd metapro_tutorial  "><pre class="notranslate"><code>cd /media/data/workspace  
mkdir metapro_tutorial/  
cd metapro_tutorial  
</code></pre></div></div>
<p dir="auto">Launch MetaPro using Singularity in interactive mode (using the command: singularity shell path/to/pipline metapro.sif</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="singularity shell /media/data/CourseData/Module6/Tools/metapro_develop.sif"><pre class="notranslate"><code>singularity shell /media/data/CourseData/Module6/Tools/metapro_develop.sif
</code></pre></div></div>
<p dir="auto"><br><br>To access MetaPro using Docker, we include instructions below: <em>(not required for the current tutorial)</em></p>
<p dir="auto">Install Docker:<br>
<a href="https://www.docker.com/products/docker-desktop" rel="nofollow">https://www.docker.com/products/docker-desktop</a></p>
<p dir="auto">Next, pull the MetaPro docker image:<br>
   docker pull parkinsonlab/metapro:develop</p>
<p dir="auto">Launch MetaPro within the Docker interactive mode:</p>
<p dir="auto">    docker run -it -v : </p>
<p dir="auto">    An example would be:<br>
    docker run -it -v /home/ubuntu/MetaPro_tutorial:/MetaPro_docker_tutorial parkinsonlab/metapro:develop</p>
<h3 dir="auto"><a id="user-content-download-the-data" class="anchor" aria-hidden="true" href="#download-the-data"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-download-the-data" href="#download-the-data"></a>ii. Download the data</h3>
<p dir="auto">Our data set consists of 150 bp single-end Illumina reads generated from mouse colon contents. Download the data and precomputed files:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="wget https://github.com/ParkinsonLab/MetaPro_tutorial/releases/download/1.0/tutorial_files.tar.gz"><pre class="notranslate"><code>wget https://github.com/ParkinsonLab/MetaPro_tutorial/releases/download/1.0/tutorial_files.tar.gz
</code></pre></div></div>
<p dir="auto">Unzip the data folder, and view contents:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="tar -xzvf tutorial_files.tar.gz  
ls"><pre class="notranslate"><code>tar -xzvf tutorial_files.tar.gz  
ls
</code></pre></div></div>
<p dir="auto">The downloaded contents include:</p>
<ul dir="auto">
<li>the input sequence file <code>mouse1.fastq</code></li>
<li>the MetaPro configuration file containing paths to tools and databases, as well as parameters <code>config_mouse_tutorial.ini</code> <em>(see below for description)</em></li>
<li>a folder containing databases required for the tutorial <code>databases\</code></li>
<li>the output directory containing some precomputed files <code>mouse1_run\</code></li>
<li>an example Cytoscape file which may be generated with MetaPro to view microbial metabolic pathway activity <code>Example.cys</code></li>
</ul>
<p dir="auto">MetaPro's tools may take a long time to run if the user does not have the necessary computing resources.  Therefore, we provide <em>pre-computed output files</em> (within the <code>mouse1_run/</code> folder) so that the user is not forced to run computationally intensive steps during the tutorial.</p>
<p dir="auto">Change the permissions of the <code>mouse1_run/</code> folder in order to view its contents through a browser ( via your IP address, http://[ public ipv4 ] ):</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="chmod -R 777 mouse1_run"><pre class="notranslate"><code>chmod -R 777 mouse1_run
</code></pre></div></div>
<p dir="auto">Access your workspace in a web browser to view the output directory: <em>(keep the tab open!)</em></p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="http://[ insert your IPv4 ]"><pre class="notranslate"><code>http://[ insert student number ].uhn-hpc.ca/
</code></pre></div></div>
<h3 dir="auto"><a id="user-content-edit-the-configuration-file" class="anchor" aria-hidden="true" href="#edit-the-configuration-file"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-edit-the-configuration-file" href="#edit-the-configuration-file"></a>iii. Edit the configuration file</h3>
<p dir="auto">MetaPro controls many of its features with a configuration file. It contains paths to tools and databases used by the program. These have already been pre-set for you in a modified config file located in the course data folder.</p>
<p dir="auto">Copy the config file and preview the contents:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="cp /media/data/CourseData/Module6/Data/config_mouse_tutorial.ini ."><pre class="notranslate"><code>cp /media/data/CourseData/Module6/Data/config_mouse_tutorial.ini .
</code></pre></div></div>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="less config_mouse_tutorial.ini"><pre class="notranslate"><code>less config_mouse_tutorial.ini
</code></pre></div></div>
<li>verify that the path to the database folder is correct.</li>

<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="[Databases]
database_path: /media/data/workspace/metapro_tutorial/databases
UniVec_Core: %(database_path)s/univec_core/UniVec_Core.fasta #(the UniVec core database)
Adapter: %(database_path)s/Trimmomatic_adapters/TruSeq3-PE-2.fa #(The adapters database)
Host: %(database_path)s/Mouse_cds/Mouse_cds.fasta #(The host database)
Rfam: %(database_path)s/Rfam/Rfam.cm #(The Infernal rRNA database)
DNA_DB: %(database_path)s/ChocoPhlAn/ChocoPhlAn.fasta #(The BWA database.  It assumes the index is in the same directory)
DNA_DB_Split: %(database_path)s/ChocoPhlAn/ChocoPhlAn_split/ #(The split database, for BLAT)
..."><pre class="notranslate"><code>[Databases]
database_path: /media/data/workspace/metapro_tutorial/databases
UniVec_Core: %(database_path)s/univec_core/UniVec_Core.fasta #(the UniVec core database)
Adapter: %(database_path)s/Trimmomatic_adapters/TruSeq3-PE-2.fa #(The adapters database)
Host: %(database_path)s/Mouse_cds/Mouse_cds.fasta #(The host database)
Rfam: %(database_path)s/Rfam/Rfam.cm #(The Infernal rRNA database)
DNA_DB: %(database_path)s/ChocoPhlAn/ChocoPhlAn.fasta #(The BWA database.  It assumes the index is in the same directory)
DNA_DB_Split: %(database_path)s/ChocoPhlAn/ChocoPhlAn_split/ #(The split database, for BLAT)
...
</code></pre></div></div>
<li>press <code>q</code> to exit</li>
<h3 dir="auto"><a id="user-content-databases-and-licenses" class="anchor" aria-hidden="true" href="#databases-and-licenses"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-databases-and-licenses" href="#databases-and-licenses"></a>iv. Databases and licenses</h3>
<p dir="auto">This tutorial relies on a few external databases and libraries to perform the filtering tasks associated with MetaPro.<br>
We have assembled the smaller databases in our precomputed files package</p>
<p dir="auto"><a href="https://ftp.ncbi.nih.gov/pub/UniVec/" rel="nofollow">The UniVec Core database</a><br>
<a href="http://ftp.ensembl.org/pub/current_fasta/mus_musculus/cds/" rel="nofollow">A mouse host sequence database</a>  In this tutorial, we will use one from Ensembl</p>
<p dir="auto">There are optional databases that are mentioned in this tutorial.  Due to their size, they will not be used here. It is highly recommended, however, that they be downloaded and used during a proper MetaPro run:</p>
<ul dir="auto">
<li><a href="http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/" rel="nofollow">The ChocoPhlan Pangenome Database</a></li>
<li>The NCBI Non-redundant (NR) Protein Database</li>
<li><a href="https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/" rel="nofollow">The GB version of the nucleotide accession2taxid table</a></li>
<li><a href="https://ccb.jhu.edu/software/centrifuge/manual.shtml#nt-database" rel="nofollow">The Centrifuge NT database</a>
<ul dir="auto">
<li>To complete the centrifuge database install, the required utilities are placed in /pipeline_tools/centrifuge</li>
</ul>
</li>
<li><a href="https://github.com/bioinformatics-centre/kaiju">The Kaiju Database</a>
<ul dir="auto">
<li>To complete the kaiju database install, the required utilities are placed in /pipeline_tools/kaiju</li>
<li>MetaPro relies on the full database. <code>(makeDB.sh -r)</code></li>
</ul>
</li>
<li><a href="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/" rel="nofollow">The NCBI Taxdump database</a></li>
<li><a href="https://www.uniprot.org/downloads" rel="nofollow">The Swiss-Prot database (fasta)</a></li>
<li><a href="http://priam.prabi.fr/REL_JAN18/Distribution.zip" rel="nofollow">The PRIAM Database</a></li>
</ul>
<p dir="auto">These optional databases require indexing prior to use.</p>
<ul dir="auto">
<li>MetaPro requires 2 versions of ChocoPhlan:
<ul dir="auto">
<li>One with all of the sequences in a single file, for BWA</li>
<li>One with all of the sequences separated, for BLAT</li>
</ul>
</li>
</ul>
<p dir="auto">The commands used to build the indexed databases are as follows <strong>[DO NOT RUN]</strong></p><ul dir="auto">
<li>To unpack chocophlan and combine all of the sequences:
<ul dir="auto">
<li>tar -xvf chocophlan.tar.gz &amp;&amp; cd chocophlan</li>
<li>for i in $(ls | grep ".gz"); do gunzip $i; done</li>
<li>for i in $(ls | grep ".m8"); do cat $i &gt;&gt; chocophlan_full.fasta</li>
<li>mv chocophlan_full.fasta ..</li>
</ul>
</li>
<li>To prepare ChocoPhlAn for BWA:  
<ul dir="auto">
<li>bwa index -a bwtsw /path/to/chocophlan_full.fasta </li>
<li>samtools faidx /path/to/chocophlan_full.fasta </li>
</ul>
</li>
<li>To prepare NR for DIAMOND:  
<ul dir="auto">
<li>diamond makedb -p 8 --in /path/to/nr -d /path/to/nr </li>
</ul>
</li>           
<li>To prepare Kaiju:
<ul dir="auto">
<li>/pipeline_tools/kaiju/makeDB.sh -r /destination/path/for/KaijuDB </li>
</ul>
</li>
<p dir="auto">Additionally, a license from MetaGeneMark is required to run to the contig assembly step</p>
<ul dir="auto">
<li><a href="http://exon.gatech.edu/Genemark/license_download.cgi" rel="nofollow">MetaGeneMark</a></li>
</ul>
<h3 dir="auto"><a id="user-content-store-paths-to-the-configuration-file-and-the-output-folder-in-variables" class="anchor" aria-hidden="true" href="#store-paths-to-the-configuration-file-and-the-output-folder-in-variables"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-store-paths-to-the-configuration-file-and-the-output-folder-in-variables" href="#store-paths-to-the-configuration-file-and-the-output-folder-in-variables"></a>v. Store paths to the configuration file and the output folder in variables</h3>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="config=/media/data/workspace/metapro_tutorial/config_mouse_tutorial.ini
output=/media/data/workspace/metapro_tutorial/mouse1_run"><pre class="notranslate"><code>config=/media/data/workspace/metapro_tutorial/config_mouse_tutorial.ini
output=/media/data/workspace/metapro_tutorial/mouse1_run
</code></pre></div></div>
<p dir="auto">These paths will be used in running MetaPro commands, structured as follows: <code>python3 /pipeline/MetaPro.py -c $config -s [sequence file] -o $output --tutorial [processing step]</code></p>
<p dir="auto">Verify the variables:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="echo $config
echo $output"><pre class="notranslate"><code>echo $config
echo $output
</code></pre></div></div>
<p dir="auto"><strong>Notes:</strong></p>
<p dir="auto">All MetaPro steps share the same file directory scheme:</p>
<ul dir="auto">
<li>data: location of interim files for each step.</li>
<li>final results: location of final outputs for each step.</li>
<li>All of MetaPro's commands are output in separate shell scripts in the relevant folders.</li>
</ul>
<h3 dir="auto"><a id="user-content-inspect-the-sequences-and-check-read-quality-with-fastqc" class="anchor" aria-hidden="true" href="#inspect-the-sequences-and-check-read-quality-with-fastqc"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-inspect-the-sequences-and-check-read-quality-with-fastqc" href="#inspect-the-sequences-and-check-read-quality-with-fastqc"></a>vi. Inspect the sequences and check read quality with FastQC</h3>
<p dir="auto">Inspect the sequences:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="less mouse1.fastq"><pre class="notranslate"><code>less mouse1.fastq
</code></pre></div></div>
<p dir="auto"><strong>Notes:</strong></p>
<ul dir="auto">
<li>Type <code>q</code> to exit <code>less</code>.</li>
</ul>
<p dir="auto">Check the read quality with FastQC:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="/pipeline_tools/FastQC/fastqc mouse1.fastq"><pre class="notranslate"><code>/pipeline_tools/FastQC/fastqc mouse1.fastq
</code></pre></div></div>
<p dir="auto">The FastQC report is generated as an HTML file <code>mouse1_fastqc.html</code>. A zip file is also generated which includes data files used to generate the report.</p>
<p dir="auto">Open the FastQC report in a browser!</p>
<p dir="auto">You can find the following information in the report:</p>
<ul dir="auto">
<li>Basic Statistics: Basic information of the mouse RNA-seq data, e.g. the total number of reads, read length, GC content.</li>
<li>Per base sequence quality: An overview of the range of quality values across all bases at each position.</li>
<li>Per Base Sequence Content: A plot showing nucleotide bias across sequence length.</li>
<li>Adapter Content: Provides information on the level of adapter contamination in your sequence sample.</li>
</ul>
<h2 dir="auto"><a id="user-content-process-the-reads" class="anchor" aria-hidden="true" href="#process-the-reads"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-process-the-reads" href="#process-the-reads"></a>Process the Reads</h2>
<h3 dir="auto"><a id="user-content-step-1-remove-adapter-sequences-and-trim-low-quality-sequences" class="anchor" aria-hidden="true" href="#step-1-remove-adapter-sequences-and-trim-low-quality-sequences"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-1-remove-adapter-sequences-and-trim-low-quality-sequences" href="#step-1-remove-adapter-sequences-and-trim-low-quality-sequences"></a>Step 1: Remove adapter sequences and trim low quality sequences.</h3>
<p dir="auto">In the first step, MetaPro removes adaptor sequences, trims low-quality reads, and removes duplicate reads in one pass.</p>
<p dir="auto">The format of the MetaPro command is:</p>
<p dir="auto">    read1='&lt;path to input sequence&gt;'<br>
    python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial quality</p>
<p dir="auto">Run the command:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="read1=/media/data/workspace/metapro_tutorial/mouse1.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial quality"><pre class="notranslate"><code>read1=/media/data/workspace/metapro_tutorial/mouse1.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial quality
</code></pre></div></div>
<p dir="auto">In this Quality-filtering stage, MetaPro will perform several actions:</p>
<ul dir="auto">
<li>filter reads below a quality score of 75</li>
<li>filter reads below a minimum length of 30 bp</li>
<li>remove adapters</li>
<li>remove duplicate reads within the dataset</li>
</ul>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 1.1: How many low quality sequences have been removed?</strong></em></p>
</blockquote>
<p dir="auto"><br><br>
Use FastQC to check the quality of the reads filtered for low quality bases and short length:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="/pipeline_tools/FastQC/fastqc mouse1_run/quality_filter/data/4_quality_filter/singletons_hq.fastq
"><pre class="notranslate"><code>/pipeline_tools/FastQC/fastqc mouse1_run/quality_filter/data/4_quality_filter/singletons_hq.fastq
</code></pre></div></div>
<p dir="auto">Navigate to the ~/4_quality_filter/ directory in your browser to view the HTML report.</p>
<p dir="auto">Compare with the previous report to see changes in the following sections:</p>
<ul dir="auto">
<li>Basic Statistics</li>
<li>Per base sequence quality</li>
</ul>
<p dir="auto"><br><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 1.2: How has the per read sequence quality curve changed in the final filtered output?</strong></em></p>
</blockquote>
<p dir="auto"><br><br>
Use FastQC to check the quality of the final filtered output:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="cd mouse1_run/quality_filter/final_results/
/pipeline_tools/FastQC/fastqc singletons_hq.fastq
cd ../../.."><pre class="notranslate"><code>cd mouse1_run/quality_filter/final_results/
/pipeline_tools/FastQC/fastqc singletons_hq.fastq
cd ../../..
</code></pre></div></div>
<p dir="auto">Compare with the previous reports to see changes in the following sections:</p>
<ul dir="auto">
<li>Basic Statistics</li>
<li>Per base sequence quality</li>
<li>Per sequence quality</li>
</ul>
  <br>
<p dir="auto"><strong>Notes on read quality filtering</strong></p>
<p dir="auto">AdapterRemoval, which was used to remove the adapters and trim low quality bases in the reads, uses a sliding window method to remove contigous regions of low quality bases in reads. However, it is worthwhile to impose an overall read quality threshold to ensure that all reads being used in our analyses are of sufficiently error-free. For this we use the tool VSEARCH which can be found at this <a href="https://github.com/torognes/vsearch">website</a> (when processing paired-end data, this step should come <strong>after</strong> the read merging step):</p>
<p dir="auto"><em>Example command only</em> <strong>[DO NOT RUN]</strong>. MetaPro already automatically perfoms this task.<br>
vsearch --fastq_filter mouse1_trim.fastq --fastq_maxee 2.0 --fastqout mouse1_qual.fastq</p>
<p dir="auto">The command line parameters are:
-   <code>--fastq_filter </code> Instructs VSEARCH to use the quality filtering algorithm to remove low quality reads
-   <code>--fastq_maxee 2.0</code> The expected error threshold. Set at 1. Any reads with quality scores suggesting that the average expected number of errors in the read are greater than 1 will be filtered.
-   <code>--fastqout</code> Indicates the output file contains the quality filtered reads</p>
 <br>
<h3 dir="auto"><a id="user-content-step-2-remove-duplicate-reads" class="anchor" aria-hidden="true" href="#step-2-remove-duplicate-reads"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-2-remove-duplicate-reads" href="#step-2-remove-duplicate-reads"></a>Step 2. Remove duplicate reads</h3>
<p dir="auto">To significantly reduce the amount of computating time required for identification and filtering of rRNA reads, we perform a dereplication step to remove duplicated reads using the software tool CD-HIT which can be obtained from this <a href="https://github.com/weizhongli/cdhit">website</a>.</p>
<p dir="auto"><em>Example command only</em> <strong>[DO NOT RUN]</strong>. MetaPro already calls this command as part of the Quality-filtering step<br>
/usr/local/prg/cd-hit-v4.6.7-2017-0501/cd-hit-auxtools/cd-hit-dup -i mouse1_qual.fastq -o mouse1_unique.fastq</p>
<p dir="auto"><strong>Notes</strong>:</p>
<p dir="auto">The command line parameters are:
-   <code>-i</code>: The input fasta or fastq file.
-   <code>-o</code>: The output file containing dereplicated sequences, where a unique representative sequence is used to represent each set of sequences with multiple replicates.
A second output file <code>mouse1_unique.fastq.clstr</code> is created which shows exactly which replicated sequences are represented by each unique sequence in the dereplicated file and a third, empty, output file, <code>mouse1_unique.fastq2.clstr</code> is also created which is only used for paired-end reads.</p>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 2.1: Can you find how many unique reads there are?</strong></em></p>
</blockquote>
<p dir="auto"><br><br>Navigate to <code>mouse1_run/quality_filter/final_results/</code> to view the FastQC report, or look at the generated <code>singletons.fastq</code> file itself in the output directory.</p>
<p dir="auto">While the number of removed duplicate reads in this small dataset is relatively low, this step can reduce file size by as much as 50-80% in larger datasets.</p>
<br>
<h3 dir="auto"><a id="user-content-step-3-remove-vector-contamination" class="anchor" aria-hidden="true" href="#step-3-remove-vector-contamination"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-3-remove-vector-contamination" href="#step-3-remove-vector-contamination"></a>Step 3. Remove vector contamination</h3>
<p dir="auto">To identify and filter reads from sources of vector, adapter, linker, and primer contamination we use the Burrows Wheeler aligner (BWA) and the BLAST-like alignment tool (BLAT) to search against a database of cow sequences. As a reference database for identifying contaminating vector and adapter sequences we rely on the UniVec_Core dataset which is a fasta file of known vectors and common sequencing adapters, linkers, and PCR Primers derived from the NCBI Univec Database.</p>
<p dir="auto">The format of the MetaPro command to perform vector filtering is:</p>
<p dir="auto">    read1='&lt;path to your quality filter final results mouse.fastq&gt;'<br>
    python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial vector</p>
<p dir="auto">Run the command:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="read1=/media/data/workspace/metapro_tutorial/mouse1_run/quality_filter/final_results/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial vector"><pre class="notranslate"><code>read1=/media/data/workspace/metapro_tutorial/mouse1_run/quality_filter/final_results/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial vector
</code></pre></div></div>
<p dir="auto">MetaPro will automatically run the following:</p>
<ul dir="auto">
<li>Index the vector contamination database.<br>
    bwa index -a bwtsw UniVec_Core<br>
    samtools faidx UniVec_Core<br>
    makeblastdb -in UniVec_Core -dbtype nucl</li>
<li>Use BWA to filter out reads aligning to the vector contamination database.<br>
    bwa mem -t 4 UniVec_Core mouse1_unique.fastq &gt; mouse1_univec_bwa.sam<br>
    samtools view -bS mouse1_univec_bwa.sam &gt; mouse1_univec_bwa.bam<br>
    samtools fastq -n -F 4 -0 mouse1_univec_bwa_contaminats.fastq mouse1_univec_bwa.bam<br>
    samtools fastq -n -f 4 -0 mouse1_univec_bwa.fastq mouse1_univec_bwa.bam</li>
<li>Use BLAT to filter out any remaining reads that align to the vector contamination database. Sicne BLAT only accepts fasta files, VSEARCH is used to first convert the reads from fastq to fasta.<br>
    vsearch --fastq_filter mouse1_univec_bwa.fastq --fastaout mouse1_univec_bwa.fasta<br>
    blat -noHead -minIdentity=90 -minScore=65  UniVec_Core mouse1_univec_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_univec.blatout</li>
<li>Lastly, a python script is used to filter the reads that BLAT does not confidently align to any vector sequences.<br>
    python3 /pipeline/Scripts/read_BLAT_Filter_v3.py single high mouse1_univec_bwa.fastq mouse1_univec.blatout mouse1_univec_blat.fastq mouse1_univec_blat_contaminats.fastq</li>
</ul>
<p dir="auto"><strong>Notes</strong>:</p>
<p dir="auto">The commands do the following tasks:<br>
-   <code>bwa index, samtools faidx, and makeblastdb</code>: Index the UniVec core database for BWA and BLAT<br>
-   <code>bwa mem</code>: Generates alignments of reads to the vector contaminant database<br>
-   <code>samtools view</code>: Converts the .sam output of bwa into .bam for the following steps<br>
-   <code>samtools fastq</code>: Generates fastq outputs for all reads that mapped to the vector contaminant database (<code>-F 4</code>) and all reads that did not map to the vector contaminant database (<code>-f 4</code>)</p>
<p dir="auto">The command line parameters for BLAT are:<br>
-   <code>-noHead</code>: Suppresses .psl header (so it's just a tab-separated file).<br>
-   <code>-minIdentity</code>: Sets minimum sequence identity is 90%.<br>
-   <code>-minScore</code>: Sets minimum score is 65. This is the matches minus the mismatches minus some sort of gap penalty.<br>
-   <code>-fine</code>: For high-quality mRNAs.<br>
-   <code>-q</code>: Query type is RNA sequence.<br>
-   <code>-t</code>: Database type is DNA sequence.</p>
<p dir="auto">The argument structure for the final python script is:<br>
    read_BLAT_Filter_v3.py &lt;operating mode: either "single" or "paired"&gt; &lt;filter stringency.  to handle paired-read conflicts.  "low" or "high"≶ &lt;Input_Reads.fq&gt; &lt;BLAT_Output_File≶ &lt;Unmapped_Reads_Output&gt; &lt;Mapped_Reads_Output&gt;`</p>
<p dir="auto">Here, BLAT does not identify any additional sequences which align to the vector contaminant database. However, we have found that BLAT is often able find alignments not identified by BWA, particularly when searching against a database consisting of whole genomes.</p>
<p dir="auto">In handling paired-ended data, cases will arise where one read maps to a vector, while the pair does not. The filter stringency decides how to resolve such cases:</p>
<ul dir="auto">
<li>Low filter stringency will only remove reads where both pairs aligned to a vector.</li>
<li>High filter stringency will remove reads where either pair aligned to a vector.</li>
</ul>
<p dir="auto"><br><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 3.1: Can you find how many reads BWA mapped to the vector database?</strong></em></p>
</blockquote>
<p dir="auto"><br></p>
<h3 dir="auto"><a id="user-content-step-4-remove-host-reads" class="anchor" aria-hidden="true" href="#step-4-remove-host-reads"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-4-remove-host-reads" href="#step-4-remove-host-reads"></a>Step 4. Remove host reads</h3>
<p dir="auto">To identify and filter host reads (here, reads of mouse origin) we repeat the steps above using a database of mouse DNA sequences. For our purposes we use a mouse genome database downloaded from Ensembl.</p>
<p dir="auto">The format of the MetaPro command is:</p>
<p dir="auto">    read1='&lt;path to your vector filter final results mouse.fastq&gt;'<br>
    python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial host</p>
<p dir="auto">Run the command:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="read1=/media/data/workspace/metapro_tutorial/mouse1_run/vector_read_filter/final_results/singletons.fastq  
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial host  "><pre class="notranslate"><code>read1=/media/data/workspace/metapro_tutorial/mouse1_run/vector_read_filter/final_results/singletons.fastq  
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial host  
</code></pre></div></div>
<p dir="auto">This call will perform the following steps:</p>
<ul dir="auto">
<li>Prepare the host database for alignment (BWA + BLAT)</li>
<li>perform alignment using BWA</li>
<li>convert the unaligned reads from BWA to a format for BLAT</li>
<li>perform alignment of the unaligned reads using BLAT</li>
<li>Run a script to remove the host reads from the input sample</li>
</ul>
<p dir="auto">The following are the commands run by the script:<br>
    bwa index -a bwtsw mouse_cds.fa<br>
    samtools faidx mouse_cds.fa<br>
    makeblastdb -in mouse_cds.fa -dbtype nucl<br>
    bwa mem -t 4 mouse_cds.fa mouse1_univec_blat.fastq &gt; mouse1_mouse_bwa.sam<br>
    samtools view -bS mouse1_mouse_bwa.sam &gt; mouse1_mouse_bwa.bam<br>
    samtools fastq -n -F 4 -0 mouse1_mouse_bwa_contaminats.fastq mouse1_mouse_bwa.bam<br>
    samtools fastq -n -f 4 -0 mouse1_mouse_bwa.fastq mouse1_mouse_bwa.bam<br>
    vsearch --fastq_filter mouse1_mouse_bwa.fastq --fastaout mouse1_mouse_bwa.fasta<br>
    blat -noHead -minIdentity=90 -minScore=65  mouse_cds.fa mouse1_mouse_bwa.fasta -fine -q=rna -t=dna -out=blast8 mouse1_mouse.blatout<br>
    ./1_BLAT_Filter.py mouse1_mouse_bwa.fastq mouse1_mouse.blatout mouse1_mouse_blat.fastq mouse1_mouse_blat_contaminats.fastq</p>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 4.1: How many reads did BWA and BLAT align to the mouse host sequence database?</strong></em></p>
</blockquote>
  <br>
<p dir="auto"><em><strong>Optional:</strong></em> In your own future analyses you can choose to complete steps 3 and 4 simultaneously by combining the vector contamination database and the host sequence database using <code>cat UniVec_Core mouse_cds.fa &gt; contaminants.fa</code>. However, doing these steps together makes it difficult to tell how much of your reads came specifically from your host organism.</p>
  <br>
<h3 dir="auto"><a id="user-content-step-5-remove-abundant-rrna-sequences--do-not-run" class="anchor" aria-hidden="true" href="#step-5-remove-abundant-rrna-sequences--do-not-run"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-5-remove-abundant-rrna-sequences--do-not-run" href="#step-5-remove-abundant-rrna-sequences--do-not-run"></a>Step 5. Remove abundant rRNA sequences *** <strong>[DO NOT RUN]</strong></h3>
<p dir="auto">rRNA genes tend to be highly expressed in all samples and must therefore be screened out to avoid lengthy downstream processing times for the assembly and annotation steps. MetaPro uses <a href="https://github.com/tseemann/barrnap">Barrnap</a> and <a href="http://infernal.janelia.org/" rel="nofollow">Infernal</a>.
You could use sequence similarity tools such as BWA or BLAST for this step, but we find Infernal, albeit slower, is more sensitive as it relies on a database of covariance models (CMs) describing rRNA sequence profiles based on the Rfam database. Due to the reliance on CMs, Infernal, can take as much as 4 hours for ~100,000 reads on a single core.  In an effort to shrink the computing time, we leverage a computing cluster's multiple cores.</p>
<p dir="auto">Here, MetaPro demonstrates the case for automation. MetaPro subdivides the input data, coordinates the concurrent processes, and collects the results into one single file after all of the scanning has been complete.</p>
<p dir="auto">MetaPro will perform the following:</p>
<ul dir="auto">
<li>subdivide the input data into user-defined chunk sizes (e.g. 1000 reads).</li>
<li>Each chunk is then run independently:
<ul dir="auto">
<li>run each chunk through Barrnap</li>
<li>using the results of Barrnap, filter the data chunk into mRNA, and leftover data for further scanning.</li>
<li>run each leftover chunk through Infernal.</li>
<li>filter the Barrnap leftover chunk using the Infernal results, to get mRNA, and "other"</li>
</ul>
</li>
<li>collect all of the data pieces (Barrnap mRNA, Infernal mRNA) into mRNA, and "other"</li>
</ul>
<p dir="auto">By running things this way, the rRNA step takes 4 minutes (as recorded with a 40-core computing node with 200 GB RAM, and an rRNA chunksize of 1000 reads), but it requires significant computing power, memory, and storage space, not available on a typical desktop PC.</p>
<p dir="auto">If you were to run this on your own, you will need the RFam database.</p>
<p dir="auto">The format of the MetaPro command is:</p>
<p dir="auto">    read1='&lt;path to your host filter final results mouse.fastq&gt;'<br>
    python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial rRNA</p>
<p dir="auto">The command would look as follows: <strong>[DO NOT RUN]</strong></p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="read1=/media/data/workspace/metapro_tutorial/mouse1_run/host_read_filter/final_results/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial rRNA"><pre class="notranslate"><code>read1=/media/data/workspace/metapro_tutorial/mouse1_run/host_read_filter/final_results/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial rRNA
</code></pre></div></div>
<p dir="auto">We have provided you with the pre-computed results:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="ls mouse1_run/rRNA_filter/final_results"><pre class="notranslate"><code>ls mouse1_run/rRNA_filter/final_results
</code></pre></div></div>
<p dir="auto">Here, we only remove about a thousand reads that map to rRNA, but in some datasets rRNA may represent up to 80% of the sequenced reads.</p>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 5.1: How many rRNA sequences were identified? How many reads are now remaining?</strong></em></p>
</blockquote>
<p dir="auto"><br></p>
<h3 dir="auto"><a id="user-content-step-6-rereplication--duplicate-repopulation" class="anchor" aria-hidden="true" href="#step-6-rereplication--duplicate-repopulation"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-6-rereplication--duplicate-repopulation" href="#step-6-rereplication--duplicate-repopulation"></a>Step 6. Rereplication / duplicate repopulation</h3>
<p dir="auto">After removing contaminants, host sequences, and rRNA, we need to replace the previously removed replicate reads back in our data set.</p>
<p dir="auto">The format of the MetaPro command is:</p>
<p dir="auto">    read1='&lt;path to your rRNA filter final results mouse.fastq&gt;'<br>
    python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial repop</p>
<p dir="auto">Run the command:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="read1=/media/data/workspace/metapro_tutorial/mouse1_run/rRNA_filter/final_results/mRNA/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial repop"><pre class="notranslate"><code>read1=/media/data/workspace/metapro_tutorial/mouse1_run/rRNA_filter/final_results/mRNA/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial repop
</code></pre></div></div>
<p dir="auto">Now that we have filtered vectors, adapters, linkers, primers, host sequences, and rRNA, check read quality of putative mRNA reads with FastQC:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="/pipeline_tools/FastQC/fastqc mouse1_run/duplicate_repopulation/final_results/singletons.fastq"><pre class="notranslate"><code>/pipeline_tools/FastQC/fastqc mouse1_run/duplicate_repopulation/final_results/singletons.fastq
</code></pre></div></div>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 6.1: How many total contaminant, host, and rRNA reads were filtered out?</strong></em></p>
</blockquote>
<p dir="auto"><br></p>
<h3 dir="auto"><a id="user-content-step-7-contig-assembly----do-not-run" class="anchor" aria-hidden="true" href="#step-7-contig-assembly----do-not-run"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-7-contig-assembly----do-not-run" href="#step-7-contig-assembly----do-not-run"></a>Step 7. Contig assembly   *** <strong>[DO NOT RUN]</strong></h3>
<p dir="auto">We have now identified the putative mRNA reads, and here we assemble the reads into contigs. Previous studies have shown that assembling reads into larger contigs significantly increases the ability to annotate them to known genes through sequence similarity searches. Here we will apply the SPAdes transcript assembly algorithm to our set of putative mRNA reads.</p>
<p dir="auto">The format of the MetaPro command is:</p>
<p dir="auto">    read1='&lt;path to your rereplicated mRNA final results mouse.fastq&gt;'<br>
    python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial contigs</p>
<p dir="auto">The command would look as follows: <strong>[DO NOT RUN]</strong></p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="read1=/media/data/workspace/metapro_tutorial/mouse1_run/duplicate_repopulation/final_results/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial contigs"><pre class="notranslate"><code>read1=/media/data/workspace/metapro_tutorial/mouse1_run/duplicate_repopulation/final_results/singletons.fastq
python3 /pipeline/MetaPro.py -c $config -s $read1 -o $output --tutorial contigs
</code></pre></div></div>
<p dir="auto">The pre-computed output is located in the following directory:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="ls mouse1_run/assemble_contigs/final_results"><pre class="notranslate"><code>ls mouse1_run/assemble_contigs/final_results
</code></pre></div></div>
<p dir="auto"><strong>Notes</strong>:
In this step, MetaPro does the following:</p>
<ul dir="auto">
<li>Assemble the reads into contigs using SPAdes</li>
<li>Use MetaGeneMark to predict the genes in these contigs</li>
<li>Use BWA to align the mRNA reads against these split-contigs to find out which reads were consumed by the process.</li>
<li>Produce a relational map of split-contig and their constituent reads.</li>
</ul>
<p dir="auto"><a href="https://cab.spbu.ru/software/spades/" rel="nofollow">SPAdes</a> assembles long contigs, but MetaPro requires that each contig only represent 1 gene. Thus the need to disassemble them, using <a href="http://exon.gatech.edu/Genemark/meta_gmhmmp.cgi" rel="nofollow">MetaGeneMark</a>  MetaGeneMark requires the user to register and obtain a free license.  Thus, we have provided the results in the precomputed files package.</p>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 7.1: How many reads were used in the assembled contigs?<br>
Hint: check how many are left in the 'singletons' (unassembled reads) file.</strong></em></p>
</blockquote>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 7.2: How many genes were predicted in the assembled contigs?<br>
Hint: try using the command <code>tail contigs.fasta</code></strong></em></p>
</blockquote>
<p dir="auto"><br></p>
<h3 dir="auto"><a id="user-content-step-8-annotate-reads-to-known-genesproteins--do-not-run" class="anchor" aria-hidden="true" href="#step-8-annotate-reads-to-known-genesproteins--do-not-run"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-8-annotate-reads-to-known-genesproteins--do-not-run" href="#step-8-annotate-reads-to-known-genesproteins--do-not-run"></a>Step 8. Annotate reads to known genes/proteins *** <strong>[DO NOT RUN]</strong></h3>
<p dir="auto">Here we will attempt to infer the specific genes our putative mRNA reads originated from. In our pipeline we rely on a tiered set of sequence similarity searches of decreasing accuracy - BWA, BLAT, and DIAMOND. While BWA provides high stringency, sequence diversity that occurs at the nucleotide level results in few matches observed for these processes. Nonetheless it is quick. To account for sequence diversity that occurs at the nucleotide level, particularly in the absence of reference microbial genomes, we use a cascaded method involving two other tools: BLAT, and DIAMOND. BLAT provides a more sensitive alignment, along with quality scores to rank the matches.  DIAMOND is used to provide more sensitive peptide-based searches, which are unaffected by sequence changes between strains.</p>
<p dir="auto">Since BWA and BLAT perform nucleotide searches, we rely on the <a href="http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/" rel="nofollow">ChocoPhlan pangenome database</a> from The Huttenhower lab, which contains reference genomes for over 10000 microbial organisms in separate .ffn files.  We create a merged copy of these sequences, and index it for BWA to use.  We leave it in its separated state for BLAT to use.</p>
<p dir="auto">For DIAMOND searches we use the Non-Redundant (NR) protein database from the NCBI.</p>
<p dir="auto">This is a computationally intensive step.  We employ our subdivision strategy here, similar to our design for rRNA removal, as seen below</p>
<ul dir="auto">
<li>The mRNA read data (contigs, remaining singletons, and remaining paired reads if applicable) are split into smaller sequence files, or 'chunks', (GA_chunksize in the configuration)</li>
<li>Each chunk is sent through BWA to be aligned against the ChocoPhlan database</li>
<li>Reads not annotated by BWA are isolated, and are sent through to BLAT.</li>
<li>Reads not annotated by BLAT are isolated, and are sent through to DIAMOND.</li>
<li>At each step, a map of annotated genes (or proteins) to constituent reads is produced.</li>
</ul> 
 <p dir="auto"></p>
<p dir="auto">In all, this leaves us with (split into chunks):</p>
<ul dir="auto">
<li>a batch of BWA-annotated gene-to-read maps</li>
<li>a batch of BLAT-annotated gene-to-read maps</li>
<li>a batch of DIAMOND-annotated protein-to-read maps</li>
<li>a batch of reads unannotated by BWA, BLAT, and DIAMOND.</li>
</ul>
<p dir="auto"></p>
<p dir="auto">These files are then merged using a custom script, and output into the GA_FINAL_MERGE folder (GA = "gene annotation"). Files include:</p>
<ul dir="auto">
<li>a gene map <code>gene_map.tsv</code></li>
<li>unannotated reads <code>GA_leftover_contigs.fasta</code> <code>GA_leftover_singletons.fasta</code></li>
<li>translated protein sequences of all genes annotated over the three annotation steps <code>all_proteins.faa</code></li>
</ul>
<p dir="auto"></p>
<p dir="auto">The format of the MetaPro command is:</p>
<p dir="auto">    read1='&lt;path to your unassembled singletons.fastq&gt;'<br>
    contig='&lt;path to your contigs.fasta&gt;'<br>
    python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial GA</p>
<p dir="auto">The command would look as follows: <strong>[DO NOT RUN]</strong>  <em>This step is heavily computationally intensive.</em></p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="read1=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial GA"><pre class="notranslate"><code>read1=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial GA
</code></pre></div></div>
<p dir="auto">We have provided pre-computed outputs in the following directory:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="ls mouse1_run/GA_FINAL_MERGE/final_results"><pre class="notranslate"><code>ls mouse1_run/GA_FINAL_MERGE/final_results
</code></pre></div></div>
<p dir="auto"><strong>Notes</strong>:</p>
<ul dir="auto">
<li>
<p dir="auto">The gene annotation step will create 4 different subdirectories: GA_BWA, GA_BLAT, GA_DIAMOND, and GA_FINAL_MERGE</p>
</li>
<li>
<p dir="auto">MetaPro will take the top-quality hits of each tool to count towards annotation:</p>
<ul dir="auto">
<li>In BWA: the CIGAR string is decoded.  Any hit with less than 90% match is rejected</li>
<li>In BLAT and DIAMOND: the only reads that pass annotation are reads with all 3 conditions satisfied:
<ul dir="auto">
<li>A sequence identity score of greater than 85</li>
<li>An alignment length of greater than 65%</li>
<li>A bitscore of over 60</li>
</ul>
</li>
</ul>
</li>
<li>
<p dir="auto">MetaPro makes a number of extra considerations to account for multi-mapped reads:</p>
<ul dir="auto">
<li>Every read scanned by BWA, BLAT, and DIAMOND is affixed with a quality score suffix taken from the aligner's report</li>
<li>In BWA: the alignment score is used.</li>
<li>In BLAT and DIAMOND: the bitscore is used.</li>
<li>Should there be a case where a read is annotated to multiple genes or proteins, each tool will use their quality scores to rank the hit
<ul dir="auto">
<li>Alignment information is iterated through, from the top of the file down to the bottom.</li>
<li>The hit with the higher score is used for all cases.</li>
<li>If the scores are tied, then the incumbent hit is used.</li>
</ul>
</li>
<li>This read-ranking is done on each chunk when the gene-to-read maps are formed.</li>
<li>There are further considerations made for paired-end annotation:
<ul dir="auto">
<li>MetaPro views paired-end data as 2 copies of the same read</li>
<li>If there are disagreements between the forward and reverse read's annotation, the quality scores are used to resolve the conflict.</li>
<li>If both reads agree on the same gene or protein, the read is counted only once.</li>
<li>This paired-read conflict resolution is performed in the GA_FINAL_MERGE step.</li>
</ul>
</li>
</ul>
</li>
<li>
<p dir="auto">Unless you are running this tutorial on a computing cluster, most systems do not have enough memory to handle indexing or searching large databases like <code>ChocoPhlan</code> (19GB) and <code>nr</code> (&gt;60GB). The descriptions in this section are purely for your information. Please use our precomputed gene, protein, and read mapping files from the tar file <code>tar -xzf tutorial_files.tar.gz</code></p>
</li>
<br>
</ul>
<h3 dir="auto"><a id="user-content-step-9-taxonomic-classification--do-not-run" class="anchor" aria-hidden="true" href="#step-9-taxonomic-classification--do-not-run"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-9-taxonomic-classification--do-not-run" href="#step-9-taxonomic-classification--do-not-run"></a>Step 9. Taxonomic Classification *** <strong>[DO NOT RUN]</strong></h3>
<p dir="auto">Now that we have putative mRNA transcripts, we can begin to infer the origins of our mRNA reads. Firstly, we will attempt to use a reference based short read classifier to infer the taxonomic orgin of our reads. Here we will use <a href="https://github.com/bioinformatics-centre/kaiju">Kaiju</a>, <a href="https://ccb.jhu.edu/software/centrifuge/manual.shtml" rel="nofollow">Centrifuge</a>, and our Gene Annotation results to generate taxonomic classifications for our reads based on a reference database.
Kaiju can classify prokaryotic reads at speeds of millions of reads per minute using the proGenomes database on a system with less than 16GB of RAM (~13GB). Using the entire NCBI nr database as a reference takes ~43GB. Similarly fast classification tools require &gt;100GB of RAM to classify reads against large databases.</p>
<p dir="auto">However, Kaiju still takes too much memory for the systems in the workshop so we have precompiled the classifications, <code>mouse1_classification.tsv</code>, in the tar file <code>tutorial_files.tar.gz</code>.
Centrifuge is a lightweight rapid microbial classification engine.  It uses methods similar to BWA and the Ferrgina-Manzini (FM) index  to make quick work of assigning taxomony.</p>
<p dir="auto">The ChocoPhlan Pangenome Database contains taxonomic information that MetaPro extracts.  Kaiju, Centrifuge, and the extracted taxa are combined using <a href="https://github.com/aametwally/WEVOTE">WEVOTE</a>.  WEVOTE is the Weighted Voting Taxonomic Identification system.  It performs consensus merging of various taxa results and reconciles the taxa identification from various sources.<br>
MetaPro uses this to settle on one confident taxon amongst Kaiju, Centrifuge, and the ChocoPhlan database choices.</p>
<p dir="auto">The format of the MetaPro command is:</p>
<p dir="auto">    read1='&lt;path to your unassembled singletons.fastq&gt;'<br>
    contig='&lt;path to your contigs.fasta&gt;'<br>
    python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial TA</p>
<p dir="auto">The command would look as follows: <strong>[DO NOT RUN]</strong>  <em>MetaPro assumes the Gene Annotation step has completed.</em></p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="read1=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial TA"><pre class="notranslate"><code>read1=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial TA
</code></pre></div></div>
<p dir="auto">We have provided pre-computed results here:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="ls mouse1_run/taxonomic_annotation/final_results"><pre class="notranslate"><code>ls mouse1_run/taxonomic_annotation/final_results
</code></pre></div></div>
<p dir="auto">We can use <a href="https://github.com/marbl/Krona/wiki">Krona</a> to generate a hierarchical multi-layered pie chart summary of the taxonomic composition of our dataset.  First, the export of MetaPro's taxonomic annotations needs to be slightly modified. Run the following commands to generate a Krona output:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="python3 /pipeline/Scripts/alter_taxa_for_krona.py mouse1_run/taxonomic_annotation/final_results/taxonomic_classifications.tsv mouse1_classification.tsv
/pipeline_tools/kaiju/kaiju2krona -t databases/nodes.dmp -n databases/names.dmp -i mouse1_classification.tsv -o mouse1_classification_Krona.txt
/pipeline_tools/KronaTools/scripts/ImportText.pl -o mouse1_classification.html mouse1_classification_Krona.txt"><pre class="notranslate"><code>python3 /pipeline/Scripts/alter_taxa_for_krona.py mouse1_run/taxonomic_annotation/final_results/taxonomic_classifications.tsv mouse1_classification.tsv
/pipeline_tools/kaiju/kaiju2krona -t databases/nodes.dmp -n databases/names.dmp -i mouse1_classification.tsv -o mouse1_classification_Krona.txt
/pipeline_tools/KronaTools/scripts/ImportText.pl -o mouse1_classification.html mouse1_classification_Krona.txt
</code></pre></div></div>
<p dir="auto">View the pie chart representation of the taxonomies detected through a web browser.</p>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 9.1: What is the most abundant family in our dataset? What is the most abundant phylum?<br>
Hint: Try decreasing the <code>Max depth</code> value on the top left of the screen and/or double-clicking on specific taxa.</strong></em></p>
</blockquote>
<p dir="auto"><br></p>
<h3 dir="auto"><a id="user-content-step-10-enzyme-function-annotation--do-not-run" class="anchor" aria-hidden="true" href="#step-10-enzyme-function-annotation--do-not-run"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-10-enzyme-function-annotation--do-not-run" href="#step-10-enzyme-function-annotation--do-not-run"></a>Step 10. Enzyme Function Annotation *** <strong>[DO NOT RUN]</strong></h3>
<p dir="auto">To help interpret our metatranscriptomic datasets from a functional perspective, we rely on mapping our data to functional networks such as metabolic pathways and maps of protein complexes. Here we predict enzyme functions (EC numbers) in our dataset, such that we may later map them to KEGG metabolic pathways.</p>
<p dir="auto">MetaPro uses three tools to produce its enzyme annotations: <a href="https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btq266" rel="nofollow">DETECT</a>, <a href="http://priam.prabi.fr/" rel="nofollow">PRIAM</a>, and DIAMOND.</p>
<p dir="auto">We use DIAMOND to identify homologs of our genes/proteins in the SWISS-PROT database that have assigned enzyme functions. This is a relatively coarse and straightforward way to annotate enzyme function by homology. We also use more robust methods for enzymatic function annotation, such as our own probability density-based enzyme function annotation tool, DETECT, and the tool PRIAM. MetaPro combines the predictions of all three tools to give two answers, a lower-confidence set of enzyme predictions <code>lq_proteins.ECs_All</code>, and a higher-confidence set of predictions <code>proteins.ECs_All</code>.</p>
<p dir="auto">PRIAM is incredibly resource-intensive and slow to run.  For the sake of time, pre-computed results have been provided.</p>
<p dir="auto">The format of the MetaPro command is:</p>
<p dir="auto">    read1='&lt;path to your unassembled singletons.fastq&gt;'<br>
    contig='&lt;path to your contigs.fasta&gt;'<br>
    python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial EC</p>
<p dir="auto">The command would look as follows: <strong>[DO NOT RUN]</strong> <em>This command assumes that MetaPro has performed the Gene annotation steps.</em></p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="read1=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial EC"><pre class="notranslate"><code>read1=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial EC
</code></pre></div></div>
<p dir="auto">The pre-computed results are provided in the following directory:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="ls mouse1_run/enzyme_annotation/final_results"><pre class="notranslate"><code>ls mouse1_run/enzyme_annotation/final_results
</code></pre></div></div>
<p dir="auto"><strong>Notes:</strong><br>
<p dir="auto"></p>
<p dir="auto">MetaPro's high-confidence and low-confidence are determined as follows:</p>
<ul dir="auto">
<li>DIAMOND: Low-confidence hits have an e-value of 1e-5 or smaller; high-confidence hits have an e-value of 1e-10 or smaller</li>
<li>PRIAM: Low-confidence hits have e-values lower than 1e-5; high-confidence hits have a probability value of 0.5 or higher</li>
<li>DETECT: No separation.</li>
</ul>
<p dir="auto"></p>
<p dir="auto">MetaPro reconciles the three sets of annotations in the following manner:</p>
<ul dir="auto">
<li>Enzymes predicted by DETECT are included, followed by annotations that agree between PRIAM and DIAMOND</li>
- In the event that multiple enzymes are annotated to the same protein:<br>
<li>Each enzyme annotation comes with a probability score.</li>
<li>MetaPro includes an enzyme co-occurence database (compiled from the ENZYME database) that contains pairs of enzymes known to exist together. Using this co-occurence database, the pipeline filters out invalid predictions. If a pair of enzymes does not exist in the database, the enzyme with the higher probability score is declared the proper annotation.</li>
<li>In cases where more than two enzymes are annotated to a protein, the top two enzymes are assigned based on the probability score.</li>
</ul>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 10.1: How many high-confidence unique enzyme functions were identified in our dataset?</strong></em></p>
</blockquote>
<p dir="auto"><br></p>
<h3 dir="auto"><a id="user-content-step-11-generate-output-files--do-not-run" class="anchor" aria-hidden="true" href="#step-11-generate-output-files--do-not-run"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-11-generate-output-files--do-not-run" href="#step-11-generate-output-files--do-not-run"></a>Step 11. Generate output files *** <strong>[DO NOT RUN]</strong></h3>
<p dir="auto">We have removed low quality bases/reads, vectors, adapters, linkers, primers, host sequences, and rRNA sequences and annotated reads to the best of our ability. We will now generate summaries of the gene counts, predicted functions and taxonomies in our microbiome.</p>
<p dir="auto">MetaPro generates many output files:</p>
<ul dir="auto">
<li>An account of read numbers during various filtering and annotation steps</li>
<li>A gene expression table of counts and RPKM values</li>
<li>RPKM values of genes in the 20-most prevalent taxa in the sample</li>
<li>A Cytoscape-compatible network file</li>
<li>An enzyme superpathway heatmap to visualize the distribution of enzymes found</li>
<li>A gene/protein-to-read map of all genes and proteins identified by MetaPro, followed by its constituent reads</li>
<li>A histogram of read quality</li>
<li>A summary of all taxa identified, followed by the number of reads associated with those taxa</li>
</ul>
<p dir="auto">The format of the MetaPro command is:</p>
<p dir="auto">    read1='&lt;path to your unassembled singletons.fastq&gt;'<br>
    contig='&lt;path to your contigs.fasta&gt;'<br>
    python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial output</p>
<p dir="auto">The command would appear as follows:  <strong>[DO NOT RUN]</strong>  <em>The command assumes that MetaPro has performed the gene, taxa, and enzyme annotations.</em></p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="read1=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial output  "><pre class="notranslate"><code>read1=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/singletons.fastq
contig=/media/data/workspace/metapro_tutorial/mouse1_run/assemble_contigs/final_results/contigs.fasta
python3 /pipeline/MetaPro.py -c $config -s $read1 --contig $contig -o $output --tutorial output  
</code></pre></div></div>
<p dir="auto">We have provided the pre-computed output files generated by MetaPro. View the files:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="ls mouse1_run/outputs/final_results/"><pre class="notranslate"><code>ls mouse1_run/outputs/final_results/
</code></pre></div></div>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 11.1: How many unique enzyme activities were predicted, at low and high confidence? View the <code>read_count.txt</code> file.</strong></em></p>
</blockquote>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 11.2: Have a look at the <code>RPKM_table.tsv</code> file. What are the most highly expressed genes?</strong></em></p>
</blockquote>
<p dir="auto"><br></p>
<blockquote>
<p dir="auto"><em><strong>Question 11.3: Which taxa appear to contribute the most metabolic activity? View the <code>enzyme_superpathway_heatmap.jpg</code>.</strong></em></p>
</blockquote>
<p dir="auto"><br></p>
<h3 dir="auto"><a id="user-content-step-12-visualize-the-results-using-a-kegg-pathway-as-a-scaffold-in-cytoscape" class="anchor" aria-hidden="true" href="#step-12-visualize-the-results-using-a-kegg-pathway-as-a-scaffold-in-cytoscape"><svg class="octicon octicon-link" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.775 3.275a.75.75 0 001.06 1.06l1.25-1.25a2 2 0 112.83 2.83l-2.5 2.5a2 2 0 01-2.83 0 .75.75 0 00-1.06 1.06 3.5 3.5 0 004.95 0l2.5-2.5a3.5 3.5 0 00-4.95-4.95l-1.25 1.25zm-4.69 9.64a2 2 0 010-2.83l2.5-2.5a2 2 0 012.83 0 .75.75 0 001.06-1.06 3.5 3.5 0 00-4.95 0l-2.5 2.5a3.5 3.5 0 004.95 4.95l1.25-1.25a.75.75 0 00-1.06-1.06l-1.25 1.25a2 2 0 01-2.83 0z"></path></svg></a><a id="user-content-step-12-visualize-the-results-using-a-kegg-pathway-as-a-scaffold-in-cytoscape" href="#step-12-visualize-the-results-using-a-kegg-pathway-as-a-scaffold-in-cytoscape"></a>Step 12. Visualize the results using a KEGG Pathway as a scaffold in Cytoscape.</h3>
<p dir="auto">This step will be performed on <strong>your workstation</strong>, using a Cytoscape file output from MetaPro.</p>
<p dir="auto">To visualize our processed microbiome dataset in the context of the carbohydrate metabolism pathways, we use the network visualization tool <strong>Cytoscape</strong> together with the <code>enhancedGraphics</code> and <code>KEGGscape</code> plugins. Some useful commands for loading in networks, node attributes and changing visual properties are provided below (there are many Cytoscape tutorials available online).  Use <a href="https://github.com/cytoscape/cytoscape/releases/3.7.2/">Cytoscape 3.7.2</a> instead of the latest version.</p>
<p dir="auto"><strong>Download the metabolic pathway</strong></p>
<p dir="auto">First, download the carbohydrate metabolism pathways from KEGG onto your workstation by pasting the following addresses into a browser:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00010.xml
https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00500.xml"><pre class="notranslate"><code>https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00010.xml
https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00500.xml
</code></pre></div></div>
<p dir="auto">Alternatively, if you are using Linux, you may use the <code>wget</code> command as follows:</p>
<div dir="auto"><div class="snippet-clipboard-content notranslate position-relative overflow-auto" data-snippet-clipboard-copy-content="wget https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00010.xml
wget https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00500.xml"><pre class="notranslate"><code>wget https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00010.xml
wget https://github.com/ParkinsonLab/Metatranscriptome-Workshop/releases/download/EC/ec00500.xml
</code></pre></div></div>
<p dir="auto">You can find other <a href="http://www.genome.jp/kegg-bin/get_htext?htext=br08901.keg" rel="nofollow">pathways on KEGG</a> which can also be imported into Cytoscape by selecting the <code>Download KGML</code> option on the top of the page for each pathway.</p>
<p dir="auto">Next, download the <code>Cytoscape_network.tsv</code> file through the browser. It is located within the <code>mouse1_run/outputs/final_results</code> folder.  This file contains expression values of enzymes (as RPKMs) of all taxa detected at &gt;1% abundance within the dataset.</p>
<p dir="auto"><strong>Install the Cytoscape plugins</strong></p>
<ul dir="auto">
<li>Select <code>Apps</code> -&gt; <code>App Manager</code></li>
<li>Search for <code>enhancedGraphics</code></li>
<li>Select <code>enhancedGraphics</code> in the middle column then click <code>Install</code> in the bottom right</li>
<li>Search for <code>KEGGScape</code></li>
<li>Select <code>KEGGScape</code> in the middle column then click <code>Install</code> in the bottom right</li>
</ul>
  <p dir="auto"></p>
<p dir="auto"><strong>Import an XML from KEGG into Cytoscape</strong></p>
<ul dir="auto">
<li>Select <code>File</code> -&gt; <code>Import</code> -&gt; <code>Network</code> -&gt; <code>File...</code></li>
<li>Select the XML file, <code>ec00010.xml</code> or <code>ec00500.xml</code> and click <code>Open</code></li>
<li>Check <code>Import pathway details from KEGG Database</code> box then select <code>OK</code></li>
</ul>
  <p dir="auto"></p>
<p dir="auto"><strong>Loading a node attribute text file (.txt)</strong> - This will map attributes (e.g. expression values) to nodes in your network which you can subsequently visualize</p>
<ul dir="auto">
<li>Select <code>File</code> -&gt; <code>Import</code> -&gt; <code>Table</code> -&gt; <code>File...</code></li>
<li>Select the <code>Cytoscape_network.tsv</code> file and click <code>Open</code></li>
<li>Change the <code>Key Column for network</code> from <code>shared name</code> to <code>KEGG_NODE_LABEL</code></li>
<li>Click OK</li>
</ul>
  <p dir="auto"></p>
<p dir="auto"><strong>Visualizing your node attributes</strong></p>
<ul dir="auto">
<li>In the left <code>Control Panel</code> select the <code>Style</code> tab</li>
<li>Check the <code>Lock node width and height</code> box</li>
<li>Click the left-most box by the <code>Size</code> panel and change the default node size to 20.0</li>
<li>Click the blank box immediately to the right of the box you clicked to change the default size, change the <code>Column</code> field to <code>RPKM</code> and the <code>Mapping Type</code> field to <code>Continuous Mapping</code></li>
<li>Click the left-most box by the <code>Image/Chart 1</code> panel, switch to the <code>Charts</code> tab, Click the doughnut ring icon, and press the <code>&gt;&gt;</code> "add all" button between the two column fields before clicking apply (make sure to remove overall RPKM from the fields that are added to the doughnut ring)</li>
<li>If you do not see the <code>Image/Chart 1</code> panel, select <code>Properties</code> -&gt; <code>Paint</code> -&gt; <code>Custom Paint 1</code> -&gt; <code>Image/Chart 1</code> from the to left corner of the control panel</li>
<li>To improve the visualization you can modify colour properties under <code>Image/Chart 1</code> -&gt; <code>Charts</code> -&gt; <code>Options</code>, or modify other properties such as Label Font Size, Label Position, Fill Color, Node location, and edge properties</li>
</ul>
  <p dir="auto"></p>
<p dir="auto"><strong>Notes:</strong></p>
<ul dir="auto">
<li>A cytoscape file with node attributes precalculated is provided for your convenience, <code>tar -xzf tutorial_files.tar.gz Example.cys</code>, feel free to open it and play with different visualizations and different layouts - compare the circular layouts with the spring embedded layouts for example. If you want to go back to the original layout then you will have to reload the file.</li>
<li>Cytoscape can be temperamental. If you don't see pie charts for the nodes, they appear as blank circles, you can show these manually. Under the 'properties' panel on the left, there is an entry labeled 'Custom Graphics 1'. Double click the empty box on the left (this is for default behavior) - this will pop up a new window with a choice of 'Images' 'Charts' and 'Gradients' - select 'Charts', choose the chart type you want (pie chart or donut for example) and select the different bacterial taxa by moving them from "Available Columns" to "Selected Columns". Finally click on 'Apply' in bottom right of window (may not be visible until you move the window).</li>
</ul>
  <p dir="auto"></p>
<p dir="auto"><strong>Visualization Questions:</strong></p>
<ul dir="auto">
<li>Which genes are most highly expressed in these two systems?</li>
<li>Which taxa are responsible for most gene expression?</li>
<li>Can you identify sub-systems (groups of interacting genes) that display anomalous taxonomic profiles?</li>
<li>Think about how you might interpret these findings; for example are certain taxa responsible for a specific set of genes that operate together to fulfill a key function?</li>
<li>Can you use the gene annotations to identify the functions of these genes through online searches?</li>
<li>Think about the implications of sequence homology searches, what may be some caveats associated with interpreting these datasets?</li>
</ul>

  
