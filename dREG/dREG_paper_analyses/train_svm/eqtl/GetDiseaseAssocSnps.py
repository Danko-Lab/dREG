import urllib2

## 0 .. 469
for i in range(0,470):
	response = urllib2.urlopen("http://regulomedb.org/GWAS/phenotype_"+str(i)+".html")
	page_source = response.read()
	
	## Get the name
	ti1 = page_source.index("H1>")
	ti2 = page_source.index("</H1")
	name = page_source[(ti1+3):ti2].replace(' ', '_')
	
	## Get indices.
	while 'chr' in page_source:
		# Get chrom.
		indx1 = page_source.index('chr')
		indx2 = page_source.index('<A HREF="rs')
		pos = page_source[indx1:indx2].replace(',', '').split(':')
		chrom = pos[0]
		chromStart = pos[1].rstrip()
		
		# Get ID.
		indx3 = page_source.index('>rs')
		indx4 = page_source.index('</A><BR')
		id = page_source[(indx3+1):indx4]
		
		# Clip page.
		last  = len(page_source)
		page_source = page_source[(indx4+1):last]
		
		# Print.
		print('\t'.join([chrom, chromStart, str(int(chromStart)+1), id, name]))

