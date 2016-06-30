library(dREG)

gdm=genomic_data_model(c(10, 25, 50, 500, 1000, 5000), c(10, 10, 30, 20, 20, 20))

#C11 bigwigs
bw_p_c11_24='/local/storage/projects/mcf7tamres/data/MCF-7_C11_GDNF_24_hr_plus.bw'
bw_m_c11_24='/local/storage/projects/mcf7tamres/data/MCF-7_C11_GDNF_24_hr_minus.bw'
bw_p_c11_1='/local/storage/projects/mcf7tamres/data/MCF-7_C11_GDNF_1_hr_plus.bw'
bw_m_c11_1='/local/storage/projects/mcf7tamres/data/MCF-7_C11_GDNF_1_hr_minus.bw'
bw_p_c11_0='/local/storage/projects/mcf7tamres/data/MCF-7_C11_GDNF_0_hr_plus.bw'
bw_m_c11_0='/local/storage/projects/mcf7tamres/data/MCF-7_C11_GDNF_0_hr_minus.bw'

#B7 bigwigs
bw_p_b7_24='/local/storage/projects/mcf7tamres/data/MCF-7_B7_GDNF_24_hr_plus.bw'
bw_m_b7_24='/local/storage/projects/mcf7tamres/data/MCF-7_B7_GDNF_24_hr_minus.bw'
bw_p_b7_1='/local/storage/projects/mcf7tamres/data/MCF-7_B7_GDNF_1_hr_plus.bw'
bw_m_b7_1='/local/storage/projects/mcf7tamres/data/MCF-7_B7_GDNF_1_hr_minus.bw'
bw_p_b7_0='/local/storage/projects/mcf7tamres/data/MCF-7_B7_GDNF_0_hr_plus.bw'
bw_m_b7_0='/local/storage/projects/mcf7tamres/data/MCF-7_B7_GDNF_0_hr_minus.bw'


getTesting<-function(bw_p, bw_m, ncores=3){
    load('/local/workdir/lauren/dtfbs_examples/chr22_tf_sites.rdata')

    test=read_genomic_data(gdm, df[,c(1,2,3)], bw_p, bw_m, ncores=ncores)

    testAll=cbind(test, df$score)

    print(dim(testAll))
    return(testAll);
}


#system.time(test2015<- getTesting(bw_p_c11_0,bw_m_c11_0, ncores=32))
#save(test2015, file='output-dreg-2015.rdata')

#system.time(test2016<- getTesting(bw_p_c11_0,bw_m_c11_0, ncores=32))
#save(test2016, file='output-dreg-2016.rdata')

#load("output-dreg-2016.rdata");
#load("output-dreg-2015.rdata");
#all(test2015==test2016);
