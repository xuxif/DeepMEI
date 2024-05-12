#HG002_10:100306854-100306954_mapClipL.sam
$file_pre=$ARGV[0];
$ran_num=$ARGV[1];
@read_ref=`samtools view ${file_pre}_mapRef.sam`;
@read_clipl=`samtools view ${file_pre}_mapClipL.sam`;
@read_clipr=`samtools view ${file_pre}_mapClipR.sam`;
@del_ref=0;
$exced_size=0;
while(@read_ref+@read_clipl+@read_clipr>95)
{
    $exced_size=1;
   if(@read_ref>@read_clipl+@read_clipr)
    {
        if($del_ref==0)
        {
            $del_ref=1;
            pop(@read_ref);
        }
        else{
            $del_ref=0;
            shift(@read_ref);
        }

    }
    elsif(@read_clipl>@read_clipr)
    {
        pop(@read_clipl);
    }
    else{
        shift(@read_clipr);
    }
}
if($exced_size==1)
{
    `cat ../head_${ran_num}.sam >${file_pre}_mapRef.sam`;
    `cat ../head_${ran_num}.sam >${file_pre}_mapClipL.sam`;
    `cat ../head_${ran_num}.sam >${file_pre}_mapClipR.sam`;
    open REF,">>${file_pre}_mapRef.sam";
    open CL,">>${file_pre}_mapClipL.sam";
    open CR,">>${file_pre}_mapClipR.sam";
    foreach $read_i (@read_ref)
    {
             print REF $read_i;
    }
    foreach $read_i (@read_clipl)
    {
             print CL $read_i;
    }
    foreach $read_i (@read_clipr)
    {
             print CR $read_i;
    }
    close(REF);
    close(CL);
    close(CR);
}
