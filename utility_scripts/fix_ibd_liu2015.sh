for trait in CD IBD UC; do
    echo ">>> Processing trait ${trait}"
    zcat EUR.${trait}.gwas_info03_filtered.assoc.gz | awk 'BEGIN{OFS="\t"}
    NR==1{
        for(i=1;i<=NF;i++){
            if($i ~ /^FRQ_A_/){ frqa_col=i; ncases=substr($i,7) }
            if($i ~ /^FRQ_U_/){ frqu_col=i; ncontrols=substr($i,7) }
            if($i == "OR"){ or_col=i }
        }
        print $0, "BETA", "FRQ", "ncases", "ncontrols"
    }
    NR>1{ print $0, log($or_col), ($frqa_col+$frqu_col)/2, ncases, ncontrols }
    ' | gzip -v > EUR.${trait}.gwas_info03_filtered.assoc.withBETA_N.txt.gz
done

