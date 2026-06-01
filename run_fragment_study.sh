#!bin/bash
start=$(date +%s) #start time
echo "start time: "`date` 
sed -i 's/\r//g' sample_path_list.txt
cat sample_path_list.txt | while read line
do

    inpath=${line}
    echo ${inpath}
    #cal endmotifs
    #python /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/endmotifs_stat_v2.2.py ${inpath} 4
    # cal LS125 sam
    #python /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/split_mis10_by125_v1.1.py ${inpath}
    #cal depth
    #bash /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/chrM_depth_v1.1.sh ${inpath}
    #cal zfsd
    #python /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/concate_depth_v4_RNAprob.py ${inpath}
    #cal zfsd feature
    #python /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/fsd_feature_cal.py ${inpath}
    # cal insertsize
    #bash /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/insertisize_find_mt.sh ${inpath}
    #concate insertsize
    python /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/insertsize_concate_v3.2.py ${inpath}
done
