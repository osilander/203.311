**[Return to the Course Home Page](../index.html)**

# TLDR Portfolio Assistance Week 4
**A/Prof Olin Silander**

[Purpose](#purpose)<br>

## Purpose

TLDR summary of how to approach to Portfolio Assessment from Week 4.

1. Find a way to get the GC content and any other useful info from the fastq.gz sequence files using `fx2tab`
2. Write the output to a file; this will automatically be in tabbed column format.
3. Import it into R; `read.table()` is one possibility with a default expectation of tabbed column format
4. Make sure you know what the columns mean. You can rename them in the text file or using `col.names=c("sequence-name", "length", "etc")` within the `read.table()` line. Note that in naming variables *and* column / row names, spaces are highly discouraged.
5. Look for correlations or patterns in the data: in chrome, [How can I look for correlations in data](https://www.google.com/search?q=How+can+I+look+for+correlations+in+data&oq=How+can+I+look+for+correlations+in+data&aqs=chrome..69i57j0i22i30j0i390i650l4.206j0j7&sourceid=chrome&ie=UTF-8)
6. Produce visualisations of your findings and explain/interpret them in a caption below the visualisaltion.