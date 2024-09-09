seqfile = RevAlign_Piwil1_cichlids.phy		* Path to the alignment file
treefile = RevAlign_Piwil1_cichlidsPiwi1314.tree	* Path to the tree file
outfile = out_Piwi1314BranchSite.txt				* Path to the output file

noisy = 3	* Display moderate amount of information on the screen
verbose = 1	* Detailed output file

seqtype = 1	* Codon data
ndata = 1	* One gene alignment
icode = 0	* Universal genetic code
cleandata = 0	* Do not remove sites with ambiguity data

model = 2	* Enable 2 or more omega for branches
NSsites = 2	* Run under model M2a
CodonFreq = 7	* Use mutation-selection model
estFreq = 0	* Use observed frequencies to calculate fitness/freq pars
clock = 0	* Assume no clock
fix_omega = 0	* Enables option to estimate omega
omega = 0.5	* Initial omega value

