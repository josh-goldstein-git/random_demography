#############################
## data preparation begins ##
#############################

## There are 380 individuals grouped into regional sub-populations.
## I select out the Greeks, whose labels begin with "gre"

## read in the raw data
x <- scan("mtdna.csv", what = character())
nchar(x)
## note the strings have different length
## because there is some header information
## and then because the labels are mixed with the sequences

## separate out the sequences
xx = x[nchar(x) > 10000] ## elements that contain both label and sequences
xx.list = strsplit(xx, split = ",") ## list of elements, split into label and sequences

## now split up this list into a vector of labels and a vector of sequences
get.first.element = function(x) {x[1]}
get.second.element = function(x) {x[2]}
labels = unlist(lapply(xx.list, get.first.element))
seqs = unlist(lapply(xx.list, get.second.element))

## Now use labels to select out the 20 Greeks
s <- grepl("^gre", labels)
my.labels = labels[s]
##  [1] "gre-1"   "gre-10"  "gre-11"  "gre-12"  "gre-15"  "gre-17"  "gre-2"
##  [8] "gre-3"   "gre-4"   "gre-5"   "gre-77"  "gre-78"  "gre-79"  "gre-80"
## [15] "gre-82"  "gre-83m" "gre-84m" "gre-85m" "gre-86"  "gre-87"
my.seqs = seqs[s]
nchar(my.seqs) ## note the sequences all have same number of bases.
##  [1] 16568 16568 16568 16568 16568 16568 16568 16568 16568 16568 16568 16568
## [13] 16568 16568 16568 16568 16568 16568 16568 16568

## now put the bases in a matrix, with each column an indiviual and
## each row a base.
my.list <- strsplit(my.seqs, "")
A <- do.call(cbind, my.list)
dim(A)
## [1] 16568    20

## coding region sequences, 576-16023 (according to "Tree construction
## and haplogroup prediction" section, but not clear if this was used for
## Intrapopulation diversity
B <- A[576:16023,]

################################
## data preparation ends here ##
################################
