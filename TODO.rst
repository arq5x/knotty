---------------
To do.
---------------
  #. Update the banners to reflect new lab.
  #. Update the Makefile for knotty (not BEDTools)
  #. Implement a custom bamToFastq tool for discordants/orphans.
  #. Implement a clean interface, e.g.

     :: $ knotty clipper -i [BAM]
        $ knotty merger -i [BAM]
        $ knotty fastq -i [BAM]
     
