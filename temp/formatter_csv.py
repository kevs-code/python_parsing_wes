#!/usr/bin/python
import re
import os, sys
#os.chdir("/home/kevin/PycharmProjects/begins/05-03-2018_mrd-myeloma-test" )

i=j=k=l=p=y=z=0
#a, b, c, d, e = ([] for ani in range(5))
a=[];b=[];c=[];d=[];e=[];g=[];h=[];m=[];r=[];t=[];q=[];u=[];#will tidy
bellpot=sys.argv[1]
with open(bellpot, 'r') as f:
    for line in f:
        if re.search("INFO=<ID=ANN", line):
            line = re.sub("^.*'Allel", "Allel", line)
            line = re.sub("' \">\s*$", ";", line)
            line = re.sub(" \| ", ";", line)
            pothole=line#occurs once now have ANN header
        if re.search("ANN=", line):
            first = re.search(r'(^.*)(CALLERS=.*)(ANN=.*\s+)(GT.*)$', line)
            base = re.sub("\s+", ",", first.group(1))#start
            annot = re.sub("^ANN=","", first.group(3))#ANNls
##############################################################################
            format, allelle1, allelle2 = re.split("\s+", first.group(4))
            mow=format
            mow = re.sub(":PGT:PID:",":",mow)
            joe=mow
            toe=mow#maybe keep format as style dependent
            joe = re.sub(":", ":t1_", joe)
            joe = re.sub("^GT:","t1_GT:",joe)
            toe = re.sub("^GT:","t2_GT:",toe)
            toe = re.sub(":",":t2_",toe)
            toe = re.sub(":t2_SA_MAP_AF:t2_SA_POST_PROB","",toe)
            #####################
            if re.search(r'GT:AD:AF:DP:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:PGT:PID:SA_MAP_AF:SA_POST_PROB',format):
                allelle2 = re.sub(":.?:.*?$","",allelle2)
                allelle1 = re.sub(":\d+\|\d+:.*?:",":",allelle1)
            #####################
            seq = (allelle1, allelle2)
            head = (joe, toe)
            s = ":"
            tupee = s.join(seq)
            oukey = s.join(head)#format header
            infohead=re.sub("=.*?;",";",first.group(2))#info header
            infohead=re.sub(";$","",infohead)
            if infohead not in e:#gives info headers and field length
                if 'DECOMP' not in infohead:#removes decomposed
                    if not re.search(";RAW_MQ;ReadPosRankSum;TLOD",infohead):
                        e.append(infohead)# not decompose = 4 now 3
            if oukey not in t:
                t.append(oukey)
                new=re.split(":",oukey)
                q.append(new)
            if format not in g:
                simples=re.split(":",format)
                print(len(simples),simples)
                g.append(format)#4 formats logged in g length vardict=80, mutect=70 +62, varscan=20
#######################still messy really
            a.append(base)#start bit content "," sep
            b.append(first.group(2))#info content ";" sep
            d.append(tupee)#format content ":" sep
            setstart=re.sub(";$","",base)
            setstart=re.sub(";",":",setstart)
            setstart=re.split(",",setstart);####start with or without AF from info
            setinfo=re.sub(";$","",first.group(2))
            parse = open('parse.csv','a')
            if re.search(r'CALLERS=vardict.*;DECOMPOSED;LEN=\d+',setinfo):
                setinfo=re.sub("$",";"*10,setinfo)
                parse.write("%s\n" % setinfo)
            if re.search(r'CALLERS=mutect2.*;DECOMPOSED;LEN=\d+',setinfo):
                setinfo=re.sub("$",";"*12,setinfo)
                parse.write("%s\n" % setinfo)
            if re.search(r'RAW_MQ=[+-]?((\d+(\.\d*)?)|\.\d+)([eE][+-]?[0-9]+)?;ReadPosRankSum',setinfo):
                setinfo=re.sub("ReadPosRankSum",";;ReadPosRankSum",setinfo)
                parse.write("%s\n" % setinfo)
            if re.search(r'ReadPosRankSum=[+-]?((\d+(\.\d*)?)|\.\d+)([eE][+-]?[0-9]+)?;TLOD=',setinfo):
                setinfo=re.sub("TLOD",";TLOD",setinfo)
            setinfo=re.split(";",setinfo);####info 4 styles 3 formats
            setformat=re.split(":",tupee);####format 4 styles 3 formats
###########################################################
            #parse = open('parse.csv','a')
            annot = re.sub("LOF=\(.*\)\s+$","",annot)#LOF for now removed save pfft trouble below
                #pfft = re.findall(r';(LOF=\(.*\))\s+$',annot)
                #parse.write("%s\n" % pfft)
                #pfft = re.search(r';(LOF=\(.*\))\s+$',annot)
                #face = re.sub("\|", ":", pfft.group(1))
                #annot = re.sub("LOF=\(.*$",face,annot)
                #pfft   #;LOF=(C2|ENSG00000166278|17|0.35),(CFB|ENSG00000243649|15|0.07),(CFB|ENSG00000244255|2|1.00)
            setannot=re.split(",",annot)

            for l in range(len(setannot)):
                setannot[l] = re.sub("\|", ",", setannot[l])
                c.append(setannot[l])
                orange=setstart+setinfo
                h.append(orange)
                m.append(setformat)
            i += 1
parse.close()
###########################################
start_header="CHROM;POS;ID;REF;ALT;QUAL;FILTER;AF;"
outty = open(bellpot+'_out_full.csv','w')
mutect2 = open(bellpot+'_mutect2_full.csv','w')
vardict = open(bellpot+'_vardict_full.csv','w')
varscan = open(bellpot+'_varscan_full.csv','w')
for j in range(len(e)):
    pot=e[j].split(';')
    if len(pot)==17:
        mutpo=';'.join(pot)
        mutpot=start_header+mutpo
        mutect2.write("%s;" % mutpot)
    if len(pot)==14:
        dictpo=';'.join(pot)
        dictpot=start_header+dictpo
        vardict.write("%s;" % dictpot)
    if len(pot)==7:
        scanpot=';'.join(pot)
        scanpot=start_header+scanpot
        varscan.write("%s;" % scanpot)
for k in range(len(q)):
    if len(q[k])==22:
        mutlis=';'.join(q[k])
        mutlist=pothole+mutlis
        mutect2.write("%s\n" % mutlist)
    if len(q[k])==12:
        scanlis=';'.join(q[k])
        scanlist=pothole+scanlis
        varscan.write("%s\n" % scanlist)
    if len(q[k])==38:
        dictlis=';'.join(q[k])
        dictlist=pothole+dictlis
        vardict.write("%s\n" % dictlist)
for p in range(len(c)):
    c[p]=re.sub(";","",c[p])
    test=c[p].split(',')
    mend =h[p]+test+m[p]
    if len(m[p]) not in u:
        u.append(len(m[p]))
    if len(m[p])==22 or len(m[p])==18:
        for item in mend:
            mutect2.write("%s;" % item)
        mutect2.write("\n")
    if len(m[p])==12:
        for item in mend:
            varscan.write("%s;" % item)
        varscan.write("\n")
    if len(m[p])==38:
        for item in mend:
            vardict.write("%s;" % item)
        vardict.write("\n")
    for item in mend:
        outty.write("%s;" % item)
    outty.write("\n")
outty.close()
mutect2.close()
vardict.close()
varscan.close()