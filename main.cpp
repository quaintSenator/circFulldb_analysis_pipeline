#include "methods.h"
const std::string inputDataDir = "E:/workplace/Bio/circRNAdata/";
const std::string outputTempDir = "./";
int main(){
    get_Mouse_chromo_micro_chromosome_distribution_matlab_data();
}
void IRES_result_preprocess(){
    std::string inputFileName = "E:/workplace/Bio/circRNAdata/IRES_source_BLASTNresult.out";
    std::ifstream inputFile(inputFileName);

    std::ofstream outputFile("IRES_afterPrePr.out");
    std::string line;
    int lineCount = 1;
    std::string currentIRES;
    std::string currentSeq;
    bool waiting4queryLine = false;
    bool waiting4subjectLine = false;
    while(std::getline(inputFile, line)){
        //汇报进度
        if(lineCount % 358 == 0){
            std::cout << "wait..."<< lineCount * 100 / 35873 << "%" << std::endl;
        }
        lineCount ++;

        if(line.size() != 0){
            //空行跳过不处理
            if(line.substr(0, 7) == "Query= "){//更换当前IRES
                currentIRES = line.substr(7, line.size() - 7);
                continue;
            }
            if(line[0] == '>'){
                currentSeq = line.substr(1, line.size() - 1);
                outputFile << "Report!" << currentIRES << "matching... " << currentSeq << std::endl;
                waiting4queryLine = true;
                waiting4subjectLine = false;
                continue;
            }
            //正常有效信息行
            if(waiting4queryLine){
                if(line[0] == 'Q'){
                    outputFile << line << std::endl;
                    waiting4subjectLine = true;
                    waiting4queryLine = false;
                }
                //丢弃行
                continue;
            }
            if(waiting4subjectLine){
                if(line[0] != 'E'){
                    outputFile << line << std::endl;
                }
                else{
                    waiting4subjectLine = false;
                }
            }
        }
    }
}
void batch_RBPmap_launcher(){
    for (int x = 5; x <= 158; x++) {
        std::string command = "perl RBPmap.pl -input allseqbatch_" + std::to_string(x)
                + " -genome mouse -db mm10 db_motifs all_ human";

        std::cout << "当前x的值: " << x << std::endl;
        std::cout << "执行命令: " << command << std::endl;

        std::system(command.c_str());
    }
}
void get_human_miRNAbinding(){
    std::string filename("E:/workplace/Bio/circRNAdata/miRNA_binding_sites_in_circRNA_sequence.human.tsv");
    // 打开文件
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << filename << std::endl;
        return;
    }
    std::string line;
    const long long totalLineCount = 23531548;
    long long currentLineCount = 0;
    std::getline(file, line);//消去第一行表头行
    currentLineCount++;
    std::ofstream output_file("human_mirna_final.json");

    output_file << "[\n";
    //行内部遍历
    while (std::getline(file, line)){
        //汇报进度
        currentLineCount++;
        if(currentLineCount % 12000 == 0){
            std::cout << "wait...completed" + std::to_string(currentLineCount * 100 / totalLineCount) << "%" << std::endl;
        }
        //这个文件的列项之间是用连续的空格分开的
        std::vector<std::string> lineStrVec(0);
        int l = 0, r = 0;

        for (; r < line.size();) //内部遍历line
        {
            //如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
            switch(line[r]){
                case '\"':
                    l = r;//l位于左侧引号
                    r++;
                    while(line[r]!='\"')r++;//r会停止在下引号位置
                    lineStrVec.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
                    r += 2;
                    l = r;
                    break;
                case '\t':
                    lineStrVec.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }
        }
        //加上尾词
        lineStrVec.push_back(line.substr(l, r - l));
        output_file << "{\n";
        output_file << "\"queryID\": \"" << lineStrVec[0] << "\",\n";
        output_file << "\"isoformID\": \"" << "chr" << lineStrVec[1] << "\",\n";
        output_file << "\"score\": " << lineStrVec[2] << ",\n";
        output_file << "\"energy\": " << lineStrVec[3] << ",\n";
        output_file << "\"queryStart\": " << lineStrVec[4] << ",\n";
        output_file << "\"queryEnd\": " << lineStrVec[5] << ",\n";
        output_file << "\"isoStart\": " << lineStrVec[6] << ",\n";
        output_file << "\"isoEnd\": " << lineStrVec[7] << ",\n";
        output_file << "\"identity\": " << lineStrVec[9] << ",\n";
        output_file << "\"similarity\": " << lineStrVec[10] << "\n";

        if(currentLineCount <= totalLineCount - 1){
            output_file << "},\n";
        }
        else{
            output_file << "}\n";
        }
    }
    output_file << "]";
    std::cout << "proc end, lineCount= " << currentLineCount << std::endl;
    // 关闭文件
    file.close();
    output_file.close();
}
void all_seq_split_into_5k(){
    const int splitSize = 3000;
    std::string inputFileName("./all_mouse_encoded.fasta");
    //std::string inputFileName("E:/workplace/cppworkplace/cmake-build-debug/all_mouse_encoded.fasta");
    // 打开文件
    std::ifstream inputFile(inputFileName);

    int batchNo = 1;
    int lineNoinBatch = 1;
    int lineNo = 1;
    std::string batchNameRoot = "allseqbatch_";
    std::string line;
    std::string outputFileName = batchNameRoot + std::to_string(batchNo);
    std::ofstream outputFile = std::ofstream(outputFileName);
    while (std::getline(inputFile, line)){
        if(line == ">seq11993"){
            std::cout << "Notice" << std::endl;
        }
        if(lineNoinBatch == splitSize){//一个batch写好了，换下一个batch
            outputFile.close();
            batchNo++;
            lineNoinBatch = 1;
            outputFileName = batchNameRoot + std::to_string(batchNo);
            outputFile = std::ofstream(outputFileName);
            std::cout << "writing into batch" << batchNo << std::endl;
        }
        if(line.size() == 0)continue;

        outputFile << line << std::endl;
        if(line[0] == '>'){
            //抬头行
        }
        else{//序列行
            lineNo++;
            lineNoinBatch++;
        }
    }
}
void get_human_TIS(){
    std::string filename("E:/workplace/Bio/circRNAdata/TIS_human.csv");
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << filename << std::endl;
        return;
    }
    std::string line;
    const int totalLineCount = 64574;
    int currentLineCount = 0;

    std::getline(file, line);//消去第一行表头行
    currentLineCount++;
    std::ofstream output_file("human_TIS_final.json");

    output_file << "[\n";
    //行内部遍历
    while (std::getline(file, line)){
        //汇报进度
        currentLineCount++;
        if(currentLineCount % 645 == 0){
            std::cout << "wait...completed" + std::to_string(currentLineCount * 100 / totalLineCount) << "%" << std::endl;
        }
        //这个文件的列项之间是用连续的空格分开的
        std::vector<std::string> lineStrVec(0);
        int l = 0, r = 0;

        for (; r < line.size();) //内部遍历line
        {
            //如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
            switch(line[r]){
                case '\"':
                    l = r;//l位于左侧引号
                    r++;
                    while(line[r]!='\"')r++;//r会停止在下引号位置
                    lineStrVec.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
                    r += 2;
                    l = r;
                    break;
                case ',':
                    lineStrVec.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }
        }
        //加上尾词
        lineStrVec.push_back(line.substr(l, r - l));
        output_file << "{\n";
        /*PK, iso_key, isoform_ID, TIS_count, TIS_on_iso_pos,
     * TIS_chr,TIS_genome_pos(5~3),TIS_strand,TIS_codon */
        if(lineStrVec[2][0] == '\"'){
            output_file << "\"isoform_ID\": " << lineStrVec[2] << ",\n";
        }
        else{
            output_file << "\"isoform_ID\": \"" << lineStrVec[2] << "\",\n";
        }
        output_file << "\"TIS_count\": " << lineStrVec[3] << ",\n";
        output_file << "\"TIS_pos\": " << lineStrVec[4] << ",\n";

        std::string interval = lineStrVec[6];
        int cutPos = 0;
        for(int i = 0; i < interval.size(); i++){
            if(interval[i] == '~'){
                cutPos = i;
                break;
            }
        }
        std::string interval_start = interval.substr(0, cutPos);
        std::string interval_end = interval.substr(cutPos + 1, interval.size() - cutPos);
        output_file << "\"interval_start\": " << interval_start << ",\n";
        output_file << "\"interval_end\": " << interval_end << ",\n";
        output_file << "\"interval_seq\": \"" << lineStrVec[8] << "\"\n";

        if(currentLineCount <= totalLineCount - 1){
            output_file << "},\n";
        }
        else{
            output_file << "}\n";
        }
    }
    output_file << "]";
    std::cout << "proc end, lineCount= " << currentLineCount << std::endl;
    // 关闭文件
    file.close();
    output_file.close();
}
void get_humanBSJallInfo(){
    //以下这段程序相对更加复杂，它被设计出来用于从humaninfo和BSJ_human里
    std::string testTissueInfo = "E:/workplace/Bio/circRNAdata/BSJ_human.csv";
    //std::string testTissueInfo = "human_BSJ_all_source_test.txt";

    std::ofstream output_file("human_allBSJs_final.json");

    // 打开文件
    std::ifstream file(testTissueInfo);
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << testTissueInfo << std::endl;
        return;
    }
    std::string line;
    long long lineCount = 0;
    //const long long totalLine = 88;
    const long long totalLine = 884637;
    std::set<std::string> stored_Otherdb_str;

    std::getline(file, line);//退掉首行
    lineCount++;
    output_file << "[\n";
    while (std::getline(file, line)) {
        //split by comma
        std::vector<std::string> splitedWords(0);
        int l = 0;
        int r = 0;
        for (; r < line.size();) //内部遍历line
        {
            //如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
            switch(line[r]){
                case '\"':
                    l = r;//l位于左侧引号
                    r++;
                    while(line[r]!='\"')r++;//r会停止在下引号位置
                    splitedWords.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
                    r += 2;
                    l = r;
                    break;
                case ',':
                    splitedWords.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }
        }
        //最后一个word被丢弃了，补
        splitedWords.push_back(line.substr(l, r - l));

        //下面是BSJ_human表的表头，
        /*
         * PK,species,BSJ_ID,gene_name,gene_id,//0-4
         * isoform_count,chr,BSJ_start,BSJ_end,strand,//5-9
         * HEK293,HEK293T,HeLa,MCF-7,SH-SY5Y,//10-14
         * SKOV3,VCaP,Adipose,Adrenal_gland,Blood,//15-19
         * Brain,Cortex,Heart,Kidney,Liver,//20-24
           Lung,Prostate,Skeletal_muscle,Smooth_muscle,Testis,//25-29
           sample_count,ciri_read,ciri_isoform_count,isoc_read,isoc_isoform_count,//30-34
           gene_name_list, tau_CIRI,tau_isoc,mouse_conservation,otherdb_Xref_count,//35-49
           Xref_DB,ORF_count,IRES_count,m6A_count,TIS_count,//40-44
           miRNA_count,RBP_count,RBP_f_count,RCS_count //45,46,47,48
        */
        //依words生成JSON
        output_file << "{\n";
        output_file << "\"species\": \"" << splitedWords[1] << "\",\n";
        output_file << "\"BSJID\": \"" << splitedWords[2] << "\",\n";

        if(splitedWords[3][0] == '\"'){
            //内部有引号，我们就不加引号
            output_file << "\"genename\": " << splitedWords[3] << ",\n";
        }
        else{
            output_file << "\"genename\": \"" << splitedWords[3] << "\",\n";
        }
        if(splitedWords[4][0] == '\"'){
            //内部有引号，我们就不加引号
            output_file << "\"geneid\": " << splitedWords[4] << ",\n";
        }
        else{
            output_file << "\"geneid\": \"" << splitedWords[4] << "\",\n";
        }

        output_file << "\"isoformcount\": " << splitedWords[5] << ",\n";
        output_file << "\"strand\": \"" << splitedWords[9] << "\",\n";
        output_file << "\"HEK293\": " << splitedWords[10] << ",\n";
        output_file << "\"HEK293T\": " << splitedWords[11] << ",\n";
        output_file << "\"HeLa\": " << splitedWords[12] << ",\n";
        output_file << "\"MCF7\": " << splitedWords[13] << ",\n";
        output_file << "\"SHSY5Y\": " << splitedWords[14] << ",\n";

        output_file << "\"SKOV3\": " << splitedWords[15] << ",\n";
        output_file << "\"VCaP\": " << splitedWords[16] << ",\n";
        output_file << "\"Adipose\": " << splitedWords[17] << ",\n";
        output_file << "\"Adrenalgland\": " << splitedWords[18] << ",\n";
        output_file << "\"Blood\": " << splitedWords[19] << ",\n";
        output_file << "\"Brain\": " << splitedWords[20] << ",\n";
        output_file << "\"Cortex\": " << splitedWords[21] << ",\n";
        output_file << "\"Heart\": " << splitedWords[22] << ",\n";
        output_file << "\"Kidney\": " << splitedWords[23] << ",\n";
        output_file << "\"Liver\": " << splitedWords[24] << ",\n";
        output_file << "\"Lung\": " << splitedWords[25] << ",\n";
        output_file << "\"Prostate\": " << splitedWords[26] << ",\n";
        output_file << "\"Skeletalmuscle\": " << splitedWords[27] << ",\n";
        output_file << "\"Smoothmuscle\": " << splitedWords[28] << ",\n";
        output_file << "\"Testis\": " << splitedWords[29] << ",\n";

        output_file << "\"samplecount\": " << splitedWords[30] << ",\n";
        output_file << "\"ciriread\": " << splitedWords[31] << ",\n";
        output_file << "\"ciriisoform_count\": " << splitedWords[32] << ",\n";
        output_file << "\"isocread\": " << splitedWords[33] << ",\n";
        output_file << "\"isocisoformcount\": " << splitedWords[34] << ",\n";
        output_file << "\"genenamelist\": \"" << splitedWords[35] << "\",\n";
        output_file << "\"tauCIRI\": \"" << splitedWords[36] << "\",\n";
        output_file << "\"tauisoc\": \"" << splitedWords[37] << "\",\n";
        output_file << "\"mouseconservation\": \"" << splitedWords[38] << "\",\n";
        output_file << "\"otherdbXref_count\": " << splitedWords[39] << ",\n";

        //output_file << "Xref_DB: \"" << splitedWords[40] << "\",\n";
        if(splitedWords[40][0] == '\"'){
            //内部有引号，我们就不加引号
            output_file << "\"XrefDB\": " << splitedWords[40] << ",\n";
        }
        else{
            output_file << "\"XrefDB\": \"" << splitedWords[40] << "\",\n";
        }
        output_file << "\"ORFcount\": " << splitedWords[41] << ",\n";
        output_file << "\"IREScount\": " << splitedWords[42] << ",\n";
        output_file << "\"m6Acount\": " << splitedWords[43] << ",\n";
        output_file << "\"TIScount\": " << splitedWords[44] << ",\n";
        output_file << "\"miRNAcount\": " << splitedWords[45] << ",\n";
        output_file << "\"RBPcount\": " << splitedWords[46] << ",\n";
        output_file << "\"RBPfcount\": " << splitedWords[47] << ",\n";
        output_file << "\"RCScount\": " << splitedWords[48] << "\n";

        //有效行其实只有totalLineCount - 1, 最后一行是total - 1，
        if(lineCount < totalLine - 1){
            output_file << "},\n"; //lineCount
        }
        else{
            output_file << "}\n";
        }

        //进度提示器
        if (lineCount % 8850 == 0) {
            std::cout << "wait..." << lineCount * 100 / totalLine << "%" << std::endl;
        }
        lineCount++;
    }
    output_file << "]\n";
    std::cout << "proc end. " << std::endl;

    // 关闭文件
    file.close();
    output_file.close();
}
void get_humanExpression(){
    //以下是对humaninfo的处理程序
    //std::string filename = "BSJ_info.txt";
    std::string filename = "E:/workplace/Bio/circRNAdata/isoform_human.csv";
    /* key,species,chr,BSJ_ID,isoform_ID,
     * exon_start,exon_end,strand,len,gene_id,
     * gene_name,algorithm,HEK293_CIRI,HEK293_isoc,HEK293T_CIRI,
     * HEK293T_isoc,HeLa_CIRI,HeLa_isoc,MCF-7_CIRI,MCF-7_isoc,
     * SH-SY5Y_CIRI,SH-SY5Y_isoc,SKOV3_CIRI,SKOV3_isoc,VCaP_CIRI,
     * VCaP_isoc, Adipose_CIRI,Adipose_isoc,Adrenal_gland_CIRI,Adrenal_gland_isoc,
     * Blood_CIRI,Blood_isoc,Brain_CIRI,Brain_isoc,Cortex_CIRI,
     * Cortex_isoc,Heart_CIRI,Heart_isoc,Kidney_CIRI,Kidney_isoc,
     * Liver_CIRI,Liver_isoc,Lung_CIRI,Lung_isoc,Prostate_CIRI,
     * Prostate_isoc,Skeletal_muscle_CIRI,Skeletal_muscle_isoc,Smooth_muscle_CIRI, Smooth_muscle_isoc,
     * Testis_CIRI, Testis_isoc, order */
    // 打开文件
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << filename << std::endl;
        return;
    }
    std::string line;
    const long long totalLineCount = 1853693;
    long long currentLineCount = 0;
    std::getline(file, line);//消去第一行表头行
    currentLineCount++;
    std::ofstream output_file("human_expression_final.json");

    output_file << "[\n";
    //行内部遍历
    while (std::getline(file, line)){
        //汇报进度
        currentLineCount++;
        if(currentLineCount % 18536 == 0){
            std::cout << "wait...completed" + std::to_string(currentLineCount * 100 / totalLineCount) << "%" << std::endl;
        }
        //这个文件的列项之间是用连续的空格分开的
        std::vector<std::string> lineStrVec(0);
        int l = 0, r = 0;

        for (; r < line.size();) //内部遍历line
        {
            //如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
            switch(line[r]){
                case '\"':
                    l = r;//l位于左侧引号
                    r++;
                    while(line[r]!='\"')r++;//r会停止在下引号位置
                    lineStrVec.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
                    r += 2;
                    l = r;
                    break;
                case ',':
                    lineStrVec.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }
        }
        //加上尾词
        lineStrVec.push_back(line.substr(l, r - l));
        output_file << "{\n";
        if(lineStrVec[4][0] == '\"'){
            output_file << "\"isoformID\": " << lineStrVec[4] << ",\n";
        }
        else{
            output_file << "\"isoformID\": \"" << lineStrVec[4] << "\",\n";
        }
        output_file << "\"order\": \"" << lineStrVec[52] << "\",\n";
        output_file << "\"readArray\": [";

        for(int i = 12; i <= 50; i++){
            int readCount = std::stoi(lineStrVec[i]);
            output_file << readCount << ", ";
        }
        //尾行没有逗号
        output_file << std::stoi(lineStrVec[51]) << "]\n";

        if(currentLineCount <= totalLineCount - 1){
            output_file << "},\n";
        }
        else{
            output_file << "}\n";
        }
    }
    output_file << "]";
    std::cout << "proc end, lineCount= " << currentLineCount << std::endl;
    // 关闭文件
    file.close();
    output_file.close();
}
void specific_expression_analysis(){
    //从获取的全亚型差异表达总表isoform_mouse_alldata.txt中提取每个亚型对每个组织差异表达矩阵
    //最终形成的差异表达矩阵的纵轴是所有亚型的ID，横轴是各个组织/细胞系，且每个组织/细胞系都有两个记录，分别对应isoc方法和CIRI方法
    std::string inputFileName = "E:/workplace/Bio/circRNAdata/isoform_mouse_alldata.txt";
    std::ifstream inputFile(inputFileName);

    // 打开文件
    if (!inputFile.is_open()) {
        std::cerr << "无法打开文件 " << inputFileName << std::endl;
        return;
    }
    const int totalLineCount = 536941;
    int lineCount = 1;
    std::string line;
    std::string headLine;

    //每个sample的总reads
    std::vector<float> totalReadVec(126, 0.0);

    //每个sample的总isoforms
    std::vector<int> countIsoVec(126, 0);

    std::vector<std::string> tableTabNames;
    std::vector<std::string> collectionNames{
        "TB","CC","HP","ST"
    };
    std::vector<float> collectionTotalReads(4,0);
    #pragma region 获取表头
    std::getline(inputFile, headLine);
    split_line_into_strvec('\t', true, tableTabNames, headLine);
    #pragma endregion

    #pragma region 第一趟-获取每一batch的read总和&最大值
    while(std::getline(inputFile, line)){
        if(lineCount % 5360 == 1){
            std::cout << "counting..." << lineCount * 100 / totalLineCount << "%\n";
        }
        lineCount++;
        std::vector<std::string> lineStrVec;
        split_line_into_strvec('\t', true, lineStrVec, line);
        for(int i = 10; i <= 125; i++){
            float currentParam = std::stof(lineStrVec[i]);
            if(currentParam > 10000){
                std::cout << "BAD" << std::endl;
            }
            totalReadVec[i] += currentParam;//加和
            if(currentParam > 0.1){
                countIsoVec[i]++;
            }
        }
    }
    #pragma endregion

    #pragma region 归类读数进四个集合
    for(int i = 10; i <= 125; i++){
       /* std::cout << i << ":"<< tableTabNames[i] << " " << countIsoVec[i] << std::endl;
        std::cout << '\t' << totalReadVec[i] << std::endl;*/
        auto collectionID2Go = decideTissueCollection(i);
        collectionTotalReads[collectionID2Go] += totalReadVec[i];
    }
    for(int i = 0; i < 4; i++){
        std::cout << collectionNames[i] << " " << collectionTotalReads[i] << std::endl;
    }
    #pragma endregion

    inputFile.clear();
    inputFile.seekg(0, std::ios::beg);
    lineCount = 1;
    std::getline(inputFile, line);//跳过头行
    std::ofstream outputFile("allmouse_exp.bed");//打开output文件

#pragma region 第二趟-计算tao
    std::cout << "2 pass - get tao..." << std::endl;

    std::getline(inputFile, line);
    while(std::getline(inputFile, line)){//遍历所有亚型
        //汇报进度
        if(lineCount % 5360 == 1){
            std::cout << "writing Into Result File..." << lineCount * 100 / totalLineCount << "%\n";
        }
        lineCount++;
        std::vector<std::string> lineStrVec;
        split_line_into_strvec('\t', true, lineStrVec, line);
        float currentIsoTAO = 0.0;
        std::vector<float> currentIsoRPMs(126, 0.0);
        std::vector<float> electedIsoRPMs(4, 0.0);
        std::vector<float> currentIsoReadCounts_in_collection(4, 0.0);
        for(int i = 10; i <= 125; i++){
            float iso_readcount_in_currentSample = std::stof(lineStrVec[i]);
            //RPM = 当前亚型在样本中的读数 * million / 样本总读数
            currentIsoRPMs[i] = iso_readcount_in_currentSample * MILLION / totalReadVec[i];//?
            currentIsoReadCounts_in_collection[decideTissueCollection(i)] += iso_readcount_in_currentSample;
        }
        //某一亚型若在若干个同集合的sample中有不同的RPM，取最大的作为该亚型在该组织的RPM
        //"TB","CC","HP","ST"
        //要用到计算完成的currentIsoRPMs，所以必须分两趟125遍历
        /*for(int i = 10; i <= 125; i++){
            auto electedRef = &electedIsoRPMs[decideTissueCollection(i)];
            //*electedRef += currentIsoRPMs[i]/collectionWeights[decideTissueCollection(i)];
            if(*electedRef < currentIsoRPMs[i]){
                *electedRef = currentIsoRPMs[i];
            }
        }*/
        for(int i = 0; i < 4; i++){
            electedIsoRPMs[i] = currentIsoReadCounts_in_collection[i] * MILLION /collectionTotalReads[i];
            if(electedIsoRPMs[i] > 10000){
                std::cout << "Bad!" << std::endl;
            }
            //currentIsoReadCounts
        }
        //
        //collectionTotalReads[i]
        float Tmax = 1.0;
        for(int i = 0 ; i < 4; i++){
            electedIsoRPMs[i] += 1.0;//规避RPM=0在对数出错，所有RPM值+1
            if(Tmax < electedIsoRPMs[i]){
                Tmax = electedIsoRPMs[i];
            }
        }
        outputFile << lineStrVec[0] << " ";
        //遍历每个集合
        float sum = 0.0;
        for(int j = 0; j < 4; j++){
            //分子求和
            float addTerm = 1 - std::log2(electedIsoRPMs[j])/std::log2(Tmax);
            sum += addTerm;
        }
        currentIsoTAO = sum / 3;
        outputFile << currentIsoTAO << " ";
        for(int i = 0; i < 4; i++){
            outputFile << currentIsoReadCounts_in_collection[i];
            outputFile << " ";
        }
        for(int i = 0; i < 4; i++){
            outputFile << (electedIsoRPMs[i] - 1);
            if(i != 3) outputFile << " ";
        }
        outputFile << std::endl;
        //lineStrVec[0]:isoform_ID
    }
#pragma endregion
}
void get_most_specific_isos(){
    std::string inputFileName = "E:/workplace/cppworkplace/cmake-build-debug/allmouse_exp.bed";
    std::ifstream inputFile(inputFileName);
    // 打开文件
    if (!inputFile.is_open()) {
        std::cerr << "无法打开文件 " << inputFileName << std::endl;
        return;
    }
    const int totalLineCount = 536941;
    int lineCount = 1;
    std::string line;

    std::vector<int> chrCounts(24, 0);
    std::vector<int> chrCounts_all(24, 0);
    std::vector<int> orgCounts(4, 0);
    std::vector<int> orgCounts_all(4, 0);
    std::vector<std::string> org_words{"TB","CC","HP","ST"};

    while(std::getline(inputFile, line)) {
        if (lineCount % 5360 == 1) {
            std::cout << "counting..." << lineCount * 100 / totalLineCount << "%\n";
        }
        lineCount++;
        std::vector<std::string> lineStrVec;
        split_line_into_strvec(' ', true, lineStrVec, line);
        float tau = std::stof(lineStrVec[1]);
        float best_RPM = 0.0;
        //2 3 4 5:readcount "TB","CC","HP","ST"
        //6 7 8 9:RPM "TB","CC","HP","ST"
        int bestAt = 0;
        for(int i = 6; i <= 9; i++){
            float currentRPM = std::stof(lineStrVec[i]);
            if(currentRPM > best_RPM){
                best_RPM = currentRPM;
                bestAt = i;
            }
            if(currentRPM != 0)orgCounts_all[i - 6]++;
        }
        int current_item_belong2chrNo = decideChrFromID(lineStrVec[0]);
        if(tau >= 0.85 && best_RPM >= 10){//汇报显著
            chrCounts[current_item_belong2chrNo]++;
            orgCounts[bestAt - 6]++;
        }
        //不论亚型是否在任何一项中显著，为染色体总亚型计数
        chrCounts_all[current_item_belong2chrNo]++;
    }
    //输出染色体-亚型/显著特异表达亚型图op
    /*for(int i = 1; i <= 22; i++){
        std::cout << "chr" << i << " " << chrCounts[i] << " " << chrCounts_all[i]<< std::endl;
    }
    std::cout << "<=====================================" << std::endl;
    for(int i = 0; i < 4; i++){
        std::cout << org_words[i] << " " << orgCounts[i] << " " << orgCounts_all[i] << std::endl;
    }*/
}
void get_RBPdata_from_all_batch(int batchNo){
    //E:\workplace\cppworkplace\cmake-build-debug
    std::string inputFileName = "E:/workplace/cppworkplace/cmake-build-debug/All_Predictions" + std::to_string(batchNo) + ".txt";
    std::ifstream inputFile(inputFileName);
    // 打开文件
    if (!inputFile.is_open()) {
        std::cerr << "无法打开文件 " << inputFileName << std::endl;
        return;
    }
    const long long totalLineCount = 1478088;
    long long lineCount = 1;
    std::string line;
    std::string currentSeq;
    std::string outputFileName = "RBPoutput_batch" + std::to_string(batchNo) + ".out";
    std::ofstream outputFile(outputFileName);

    while(std::getline(inputFile, line)){
        if(lineCount % 14780 == 1){
            std::cout << "processing batch" << batchNo << " " << lineCount * 100 / totalLineCount << "%\n";
        }
        lineCount++;

        if(line.substr(0, 3) == "seq"){ //更新current seq
            currentSeq = line;
            continue;
        }
        if(line.substr(0, 6) == "Report"){
            //组装currentSeq + Report
            int l = 0, r = 0, numberPos = 0, numberEnd = 0;
            for(int i = 0; i < line.size(); i++){
                if(line[i] == ':'){
                    l = i + 2;//词头
                }
                if(line[i] == '('){
                    r = i - 1;//词尾
                }
                if(line[i] == 'n' && i + 1 < line.size() && line[i + 1] == 'g'){
                    numberPos = i + 3;
                }
                if(line[i] == 's' && i + 1 < line.size() && line[i + 1] == 'i'){
                    numberEnd = i - 1;//注意，numberEnd是一个空格
                }
            }
            std::string ProteinName = line.substr(l, r - l + 1);
            std::string ProteinCountStr = line.substr(numberPos, numberEnd - numberPos);
            outputFile << currentSeq << " " << ProteinName << " " << ProteinCountStr << std::endl;
        }
    }
}
int decideTissueCollection(int index){
    //汇总得到四个collection的总reads
    //"TB","CC","HP","ST"
    //ST:108-110, 123-125
    //HP:105-107, 120-122
    //CC:102-104, 117-119
    //TB:else
    if(index < 102)return 0;//TB
    if(index <= 104)return 1;//CC
    if(index <= 107)return 2;//HP
    if(index <= 110)return 3;//ST
    if(index <= 116)return 0;//TB
    if(index <= 119)return 1;//CC
    if(index <= 122)return 2;//HP
    return 3;//ST
}
int decideChrFromID(std::string &s){
    if(s.size() < 1)return 0;
    int cutPos = 0;
    for(int i = 0; i < s.size(); i++){
        if(s[i] == ':'){
            cutPos = i;
            break;
        }
    }
    std::string result_str = s.substr(3, cutPos - 3);
    if(result_str[0] == 'X' || result_str[0] == 'x')return 20;
    if(result_str[0] == 'Y' || result_str[0] == 'y')return 21;
    if(result_str[0] == 'M' || result_str[0] == 'm')return 22;
    if(result_str[0] == 'U' || result_str[0] == 'u')return 0;
    float result = 0;
    try{
         result = std::stoi(result_str);
    }
    catch(std::invalid_argument){
        std::cout << "Bad at " << result_str << std::endl;
    }
    return result;
}
void judge_Human_Mouse_BSJs_has_collapse(){
    //意在判断Human和Mouse是否有重复亚型
    std::string human_iso_file_name = "E:/workplace/Bio/circRNAdata/isoform_human.csv";
    std::string mouse_iso_file_name = "E:/workplace/Bio/circRNAdata/isoform_mouse_alldata.txt";

    std::ifstream human_iso_file(human_iso_file_name);
    std::ifstream mouse_iso_file(mouse_iso_file_name);

    std::set<std::string> human_iso_map;
    //human_isoform 的前几项
    /*0,Human,chr10,chr10|100036449|100036604|+,
     * chr10|100036449|100036604|+,100036449,100036604,+,
     * 156,Unknown,Unknown,CIRI-long
     * */

#pragma region 1pass-建立set
    int lineCount = 1;
    int totalLineCount = 1848000;
    std::string line;
    std::getline(human_iso_file, line);//退掉首行
    while(std::getline(human_iso_file, line)){
        if (lineCount % 18480 == 1) {
            std::cout << "counting..." << lineCount * 100 / totalLineCount << "%\n";
        }
        lineCount++;
        std::vector<std::string> lineStrVec;
        int l = 0, r = 0;
        for (; r < line.size();) //内部遍历line
        {
            switch (line[r]) {
                case '\"'://如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
                    l = r;//l位于左侧引号
                    r++;
                    while (line[r] != '\"')r++;//r会停止在下引号位置
                    lineStrVec.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
                    r += 2;
                    l = r;
                    break;
                case ',':
                    lineStrVec.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }
        }
        lineStrVec.push_back(line.substr(l, line.size() - l));
        //line read out
        std::string isoID = "";
        //[5] "\"100042282,100048758,100054347,100057013,100063614,100065188\""
        if(lineStrVec[4][0] == '\"'){
            isoID = lineStrVec[4].substr(1, lineStrVec[4].size() - 2);
        }
        else{
            isoID = lineStrVec[4];
        }
        human_iso_map.insert(isoID);
    }
#pragma endregion

#pragma region 2pass-遍历mouse_iso
    lineCount = 1;
    int collapseCount = 0;
    totalLineCount = 535920;
    std::getline(mouse_iso_file, line);//退掉头行
    //一个重大问题：isoID格式差异
    while(std::getline(mouse_iso_file, line)){
        if (lineCount % 5359 == 1) {
            std::cout << "mapping..." << lineCount * 100 / totalLineCount << "%\n";
        }
        lineCount++;
        std::vector<std::string> lineStrVec;
        std::vector<std::string> seg_begins;
        std::vector<std::string> seg_ends;
        int l = 0, r = 0;
        for (; r < line.size();) //内部遍历line
        {
            switch (line[r]) {
                case '\"'://如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
                    l = r;//l位于左侧引号
                    r++;
                    while (line[r] != '\"')r++;//r会停止在下引号位置
                    lineStrVec.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
                    r += 2;
                    l = r;
                    break;
                case '\t':
                    lineStrVec.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }
        }
        lineStrVec.push_back(line.substr(l, line.size() - l));
        //line cut read out
        std::stringstream isoID;
        isoID << lineStrVec[3] << "|";
        //从mouse数据里提取段信息，重新拼接成

        //human = "chr10|100054374,100057013,100063614|100054446,100057152,100063682|-"
        //mouse = "chr10|100054374-100054446,100057013-100057152,1000636140-100063682|+"
        //[4]"142303943,142304297,142312109,142314406,142352734"
        //[5]"142303979,142304469,142312214,142314448,142352885"
        //cut [4] & [5]
        int ll = 0;
        for(int i = 0; i < lineStrVec[4].size(); i++){
            if(lineStrVec[4][i] == ','){
                seg_begins.push_back(substr_by_be(lineStrVec[4], ll, i - 1));
                ll = i + 1;
            }
        }
        seg_begins.push_back(substr_by_be(lineStrVec[4],ll, lineStrVec[4].size() - 1));
        ll = 0;
        for(int i = 0; i < lineStrVec[5].size(); i++){
            if(lineStrVec[5][i] == ','){
                seg_ends.push_back(substr_by_be(lineStrVec[5], ll, i - 1));
                ll = i + 1;
            }
        }
        seg_ends.push_back(substr_by_be(lineStrVec[5],ll, lineStrVec[5].size() - 1));
        for(int i = 0; i < seg_begins.size() - 1; i++){
            //目标 "chr10|100054374,100057013,100063614|100054446,100057152,100063682|-"
            isoID << seg_begins[i] << ",";
        }
        isoID << seg_begins[seg_begins.size() - 1] << "|" ;//
        for(int i = 0; i < seg_ends.size() - 1; i++){
            //目标 "chr10|100054374,100057013,100063614|100054446,100057152,100063682|-"
            isoID << seg_ends[i] << ",";
        }
        isoID << seg_begins[seg_ends.size() - 1] << "|" ;//

        isoID << lineStrVec[6];
        std::string isoIDstr = isoID.str();
        if(human_iso_map.find(isoIDstr) != human_iso_map.end()){
            std::cout << isoIDstr << std::endl;
            collapseCount++;
        }
    }

    std::cout << collapseCount << std::endl;
#pragma endregion
}
void count_Human_BSJs(){
    std::string inputFileName = "E:/workplace/cppworkplace/cmake-build-debug/mouse_isos2Append.json";
    std::ifstream inputFile(inputFileName);
    std::string line;
    std::set<std::string> bsj_set;
    std::getline(inputFile, line);
    while(std::getline(inputFile, line)){
        if(line.substr(0, 7) == "\"BSJID\""){
            std::string cutBSJID = line.substr(10, line.size() - 10);
            if(bsj_set.find(cutBSJID) == bsj_set.end()){
                bsj_set.insert(cutBSJID);
            }
        }
    }
    std::cout << bsj_set.size() << std::endl;
}
void count_Mouse_BSJs(){
    std::string inputFileName = "E:/workplace/cppworkplace/cmake-build-debug/human_allBSJs_final.json";

    std::ifstream inputFile(inputFileName);
    std::string line;
    std::set<std::string> bsj_set;
    std::getline(inputFile, line);
    while(std::getline(inputFile, line)){
        if(line.substr(0, 7) == "\"BSJID\""){
            //"BSJID": "chr15:41661243-41662535",
            std::string cutBSJID = line.substr(9, line.size() - 9);
            if(bsj_set.find(cutBSJID) == bsj_set.end()){
                bsj_set.insert(cutBSJID);
            }
        }
    }
    std::cout << bsj_set.size() << std::endl;
}
void mouse_id_recoder(std::string& codingStr){
    std::vector<int> cutPosList;
    std::vector<int> midPosList;
    int segmentBegin = 0;
    for(int i = 0; i < codingStr.size(); i++){
        if(codingStr[i] == '-'){
            midPosList.push_back(i);
        }
        if(codingStr[i] == ','){
            cutPosList.push_back(i);
        }
        if(codingStr[i] == '|'){
           segmentBegin = i;
        }
    }
    //chr10| ==> segmentBegin = 5, substr_len = 6
    std::string codedResult = codingStr.substr(0, segmentBegin + 1);//|字符已经加入
}
void get_humanIso(){
    //以下是对humaninfo的处理程序
    //std::string filename = "BSJ_info.txt";
    std::string filename = "E:/workplace/Bio/circRNAdata/isoform_human.csv";
    /* key,species,chr,BSJ_ID,isoform_ID,
     * exon_start,exon_end,strand,len,gene_id,
     * gene_name,algorithm,HEK293_CIRI,HEK293_isoc,HEK293T_CIRI,
     * HEK293T_isoc,HeLa_CIRI,HeLa_isoc,MCF-7_CIRI,MCF-7_isoc,
     * SH-SY5Y_CIRI,SH-SY5Y_isoc,SKOV3_CIRI,SKOV3_isoc,VCaP_CIRI,
     * VCaP_isoc, Adipose_CIRI,Adipose_isoc,Adrenal_gland_CIRI,Adrenal_gland_isoc,
     * Blood_CIRI,Blood_isoc,Brain_CIRI,Brain_isoc,Cortex_CIRI,
     * Cortex_isoc,Heart_CIRI,Heart_isoc,Kidney_CIRI,Kidney_isoc,
     * Liver_CIRI,Liver_isoc,Lung_CIRI,Lung_isoc,Prostate_CIRI,
     * Prostate_isoc,Skeletal_muscle_CIRI,Skeletal_muscle_isoc,Smooth_muscle_CIRI, Smooth_muscle_isoc,
     * Testis_CIRI, Testis_isoc, order */
    // 打开文件
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << filename << std::endl;
        return;
    }
    std::string line;
    const long long totalLineCount = 1853693;
    long long currentLineCount = 0;
    std::getline(file, line);//消去第一行表头行
    currentLineCount++;
    std::ofstream output_file("human_iso_final.json");

    std::vector<std::string> detectionDictionary{"HEK293_CIRI","HEK293_isoc","HEK293T_CIRI",
                                                 "HEK293T_isoc","HeLa_CIRI","HeLa_isoc","MCF-7_CIRI","MCF-7_isoc",
                                                 "SH-SY5Y_CIRI","SH-SY5Y_isoc","SKOV3_CIRI","SKOV3_isoc","VCaP_CIRI",
                                                 "VCaP_isoc", "Adipose_CIRI","Adipose_isoc","Adrenal_gland_CIRI","Adrenal_gland_isoc",
                                                 "Blood_CIRI","Blood_isoc","Brain_CIRI","Brain_isoc","Cortex_CIRI",
                                                 "Cortex_isoc","Heart_CIRI","Heart_isoc","Kidney_CIRI","Kidney_isoc",
                                                 "Liver_CIRI","Liver_isoc","Lung_CIRI","Lung_isoc","Prostate_CIRI",
                                                 "Prostate_isoc","Skeletal_muscle_CIRI","Skeletal_muscle_isoc",
                                                 "Smooth_muscle_CIRI", "Smooth_muscle_isoc",
                                                 "Testis_CIRI", "Testis_isoc"};
    output_file << "[\n";
    //行内部遍历
    while (std::getline(file, line)){
        //汇报进度
        currentLineCount++;
        if(currentLineCount % 18536 == 0){
            std::cout << "wait...completed" + std::to_string(currentLineCount * 100 / totalLineCount) << "%" << std::endl;
        }
        //这个文件的列项之间是用连续的空格分开的
        std::vector<std::string> lineStrVec(0);
        int l = 0, r = 0;

        for (; r < line.size();) //内部遍历line
        {
            //如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
            switch(line[r]){
                case '\"':
                    l = r;//l位于左侧引号
                    r++;
                    while(line[r]!='\"')r++;//r会停止在下引号位置
                    lineStrVec.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
                    r += 2;
                    l = r;
                    break;
                case ',':
                    lineStrVec.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }
        }
        lineStrVec.push_back(line.substr(l, r - l));
        output_file << "{\n";
        output_file << "\"BSJID\": \"" << lineStrVec[3] << "\",\n";
        output_file << "\"species\": \"" << lineStrVec[1] << "\",\n";
        if(lineStrVec[4][0] == '\"'){
            output_file << "\"isoformID\": " << lineStrVec[4] << ",\n";
        }
        else{
            output_file << "\"isoformID\": \"" << lineStrVec[4] << "\",\n";
        }
        if(lineStrVec[5][0] == '\"'){
            output_file << "\"exonstart\": " << lineStrVec[5] << ",\n";
        }
        else{
            output_file << "\"exonstart\": \"" << lineStrVec[5] << "\",\n";
        }
        if(lineStrVec[6][0] == '\"'){
            output_file << "\"exonend\": " << lineStrVec[6] << ",\n";
        }
        else{
            output_file << "\"exonend\": \"" << lineStrVec[6] << "\",\n";
        }
        if(lineStrVec[11][0] == '\"'){
            output_file << "\"algo\": " << lineStrVec[11] << ",\n";
        }
        else{
            output_file << "\"algo\": \"" << lineStrVec[11] << "\",\n";
        }
        output_file << "\"order\": \"" << lineStrVec[52] << "\",\n";
        std::string detection = "";
        for(int i = 12; i <= 51; i++){
            //i = 12 -> vector[0]
            if(lineStrVec[i] != "0.0"){
                detection += detectionDictionary[i - 12];
            }
        }
        output_file << "\"detected\": \"" << detection << "\"\n";


        if(currentLineCount <= totalLineCount - 1){
            output_file << "},\n";
        }
        else{
            output_file << "}\n";
        }
    }
    output_file << "]";
    std::cout << "proc end, lineCount= " << currentLineCount << std::endl;
    // 关闭文件
    file.close();
    output_file.close();
}
void merge_Human_Mouse_isos(){
    //意在生成mouse_isos向human_isos的插入json，从而对human_isos直接插入形成最终总iso表
    std::string mouse_iso_file_name = "E:/workplace/Bio/circRNAdata/isoform_mouse_alldata.txt";
    std::string outputFileName = "mouse_isos2Append.json";
    std::ifstream mouse_iso_file(mouse_iso_file_name);
    std::ofstream outputJsonFile(outputFileName);
    std::string line;

    int lineCount = 1;
    int collapseCount = 0;
    int totalLineCount = 535920;
    std::getline(mouse_iso_file, line);//退掉头行
    outputJsonFile << "[\n";
    while(std::getline(mouse_iso_file, line)){
        if (lineCount % 5359 == 1) {
            std::cout << "mapping..." << lineCount * 100 / totalLineCount << "%\n";
        }
        lineCount++;
        std::vector<std::string> lineStrVec;
        std::vector<std::string> seg_begins;
        std::vector<std::string> seg_ends;
        split_line_into_strvec('\t',true, lineStrVec, line);


        /*#pragma region 重组isoID
        std::stringstream isoID;
        isoID << lineStrVec[3] << "|";

        int ll = 0;
        for(int i = 0; i < lineStrVec[4].size(); i++){
            if(lineStrVec[4][i] == ','){
                seg_begins.push_back(substr_by_be(lineStrVec[4], ll, i - 1));
                ll = i + 1;
            }
        }
        seg_begins.push_back(substr_by_be(lineStrVec[4],ll, lineStrVec[4].size() - 1));
        ll = 0;
        for(int i = 0; i < lineStrVec[5].size(); i++){
            if(lineStrVec[5][i] == ','){
                seg_ends.push_back(substr_by_be(lineStrVec[5], ll, i - 1));
                ll = i + 1;
            }
        }
        seg_ends.push_back(substr_by_be(lineStrVec[5],ll, lineStrVec[5].size() - 1));
        for(int i = 0; i < seg_begins.size() - 1; i++){
            //目标 "chr10|100054374,100057013,100063614|100054446,100057152,100063682|-"
            isoID << seg_begins[i] << ",";
        }
        isoID << seg_begins[seg_begins.size() - 1] << "|" ;//
        for(int i = 0; i < seg_ends.size() - 1; i++){
            //目标 "chr10|100054374,100057013,100063614|100054446,100057152,100063682|-"
            isoID << seg_ends[i] << ",";
        }
        isoID << seg_ends[seg_ends.size() - 1] << "|" ;//

        isoID << lineStrVec[6];
        std::string isoIDstr = isoID.str();
        #pragma endregion*/
        outputJsonFile << "{\n";
        std::string BSJID = "";
        if(lineStrVec[2][0] == '\"'){
            BSJID = lineStrVec[2].substr(1, lineStrVec[2].size() - 2);
        }
        else BSJID = lineStrVec[2];
        //小鼠BSJID的chr10:与人类BSJID的chr10|风格不符，采用人类版风格
        for(int i = 0; i < BSJID.size(); i++){
            if(BSJID[i] == ':'){
                BSJID[i] = '|';
                break;
            }
        }

        outputJsonFile << "\"BSJID\": \"" << BSJID << "\",\n";
        outputJsonFile << "\"species\": \"Mouse\",\n";
        outputJsonFile << "\"isoformID\":\"" << lineStrVec[0] << "\",\n";
        outputJsonFile << "\"exonstart\": \"" << lineStrVec[4] << "\",\n";
        outputJsonFile << "\"exonend\": \"" << lineStrVec[5] << "\"\n";

        if(lineCount == 536941){
            outputJsonFile << "}\n";
        }
        else{
            outputJsonFile << "},\n";
        }
    }
    outputJsonFile << "]\n";
    std::cout << lineCount << std::endl;
}
void merge_Human_Mouse_exps(){
    std::string human_exps_filename = inputDataDir + "isoform_human.csv";
    std::ifstream human_exps_file(human_exps_filename);

    /*
     * 当前的一个human_exps记录：
  isoformID: 'chr10|1000677,1000948,1007073|1000868,1000977,1007128|+',
  order: 'iso#002',
  readArray: [
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  ]
    10个样本，分别用两种不同的鉴定工具，共计20项；
     但是现在引入了小鼠组，小鼠的表达列加工后已经变成tau
}
     * */
    std::string line;
    //退首行
    std::getline(human_exps_file, line);
    while(std::getline(human_exps_file, line)){
        /*key,species,chr,BSJ_ID,isoform_ID,
         * exon_start,exon_end,strand,len,gene_id,
         * gene_name,algorithm,HEK293_CIRI,HEK293_isoc,HEK293T_CIRI,
         * HEK293T_isoc,HeLa_CIRI,HeLa_isoc,MCF-7_CIRI,MCF-7_isoc,
         * SH-SY5Y_CIRI,SH-SY5Y_isoc,SKOV3_CIRI,SKOV3_isoc,VCaP_CIRI,
         * VCaP_isoc,Adipose_CIRI,Adipose_isoc,Adrenal_gland_CIRI,Adrenal_gland_isoc,
         * Blood_CIRI,Blood_isoc,Brain_CIRI,Brain_isoc,Cortex_CIRI,
         * Cortex_isoc,Heart_CIRI,Heart_isoc,Kidney_CIRI,Kidney_isoc,
         * Liver_CIRI,Liver_isoc,Lung_CIRI,Lung_isoc,Prostate_CIRI,
         * Prostate_isoc,Skeletal_muscle_CIRI,Skeletal_muscle_isoc,Smooth_muscle_CIRI,Smooth_muscle_isoc,
         * Testis_CIRI,Testis_isoc,order
         * */
        std::vector<std::string> lineStrVec;
        split_line_into_strvec(',', true, lineStrVec, line);
    }
}

void get_ORFs(){
    std::string ORF_input_filename = inputDataDir + "isoform_mouse.repeat4times.fa.ORF.out";
    std::ifstream ORF_source_file(ORF_input_filename);
    std::string iso_input_filename = inputDataDir + "isoform_mouse_alldata.txt";
    std::ifstream iso_source_file(iso_input_filename);

    std::string outputFileName = "ORF_processed.bed";
    std::ofstream outputFile(outputFileName);

    std::string line;
    std::map<std::string, isoform_info> isoID_map;

    #pragma region 第一趟-建isoID_map
    //chr10:100051846-100064146|100051846-100051959,100064084-100064146 输入格式

    int lineCount = 0;
    int totalLineCount = 536800;
    std::getline(iso_source_file, line);//退掉第一行
    while(std::getline(iso_source_file, line)){
        lineCount++;
        if(lineCount % 8800 == 1){
            std::cout << "building isoID_map..." << lineCount * 100 / totalLineCount << "%" << std::endl;
        }
        std::vector<std::string> lineStrVec;
        split_line_into_strvec('\t',false, lineStrVec, line);
        isoID_map.insert(std::pair<std::string, isoform_info>(lineStrVec[0],
                                                     isoform_info(std::stoi(lineStrVec[7]), lineStrVec[6][0])));
    }
    #pragma endregion
    long long total_ORF_count = 0;
    long long hit_ORF_count = 0;
    #pragma region 第二趟-命中并填入Json
    std::string processingIsoID;
    totalLineCount = 20965208;
    lineCount = 0;
    int overlong_ORF_count = 0;
    int cORF_count = 0;
    int mem_lbit;
    int mem_rbit;
    isoform_info mem_circInfo(0,'+');
    std::string proteinName;
    bool hasHit = false;
    outputFile << "[\n";
    while(std::getline(ORF_source_file, line)){
        lineCount ++;
        if(lineCount % 209652 == 1){
            std::cout << "hitting and writing..." << lineCount * 100 / totalLineCount << "%" << std::endl;
        }
        std::string lBitPosStr;
        std::string rBitPosStr;
        //lcl|ORF26_chr10:100051846-100089271|100051846-100051959,100055765-100055858,100064084-100064146,
        // 100076735-100076905,100079974-100080130,
        // 100080857-100080940,100087347-100087456,100088153-100088220,100089196-100089271:2013:1891 unnamed protein product
        if(line[0] == '>'){//抬头行
            int ID_substr_start = 0;
            int ID_substr_end = 0;
            int lr_split = 0;
            int lr_end = 0;
            for(int i = 0; i < line.size(); i++){
                if(line[i] == '_'){
                    ID_substr_start = i + 1;
                }
            }
            //lcl|ORF26_chr10:100051846-100089271|100051846-100051959,100055765-100055858,100064084-100064146,
            // 100076735-100076905,100079974-100080130,
            // 100080857-100080940,100087347-100087456,100088153-100088220,100089196-100089271:2013:1891 unnamed protein product
            for(int i = ID_substr_start + 6; i < line.size(); i++){
                if(line[i] == ':'){//:2013:1891 当中第一个冒号
                    ID_substr_end = i - 1;
                    break;
                }
            }
            for(int i = ID_substr_end + 2; i < line.size(); i++){
                if(line[i] == ':'){
                    lr_split = i;
                    break;
                }
            }
            for(int i = lr_split + 1; i < line.size(); i++){
                if(line[i] == ' '){
                    lr_end = i;
                    break;
                }
            }
            proteinName = substr_by_be(line, lr_end + 1, line.size() - 1);
            lBitPosStr = substr_by_be(line, ID_substr_end + 2, lr_split - 1);
            rBitPosStr = substr_by_be(line, lr_split + 1, lr_end - 1);
            mem_lbit = std::stoi(lBitPosStr);
            mem_rbit = std::stoi(rBitPosStr);
            if(mem_lbit > mem_rbit)std::swap(mem_lbit, mem_rbit);
            processingIsoID = substr_by_be(line, ID_substr_start, ID_substr_end);
            //std::string prIsoIDtest1 = processingIsoID + "|+";
            //std::string prIsoIDtest2 = processingIsoID + "|-";
            auto mapFindResult = isoID_map.find(processingIsoID);
            if(mapFindResult != isoID_map.end()){
                hit_ORF_count++;
                hasHit = true;
                mem_circInfo = isoform_info(mapFindResult -> second);
            }
        }
        else{
            //蛋白质组成行
            if(!hasHit){
                continue;
            }
            else{
                bool shouldWrite = false;
                //确定命中开始判断ORF性质
                //恢复hasHit置位
                int ORFtype = judgeORFvalid(mem_circInfo.iso_len, mem_lbit, mem_rbit);
                switch(ORFtype){
                    case 0:
                        overlong_ORF_count++;
                        break;
                    case 1:
                        shouldWrite = true;
                        break;
                    case 2:
                        shouldWrite = true;
                        cORF_count++;
                        break;
                    default:
                        break;
                }
                if(shouldWrite){
                    outputFile << "{\n" << "\"isoformID\": \"" << processingIsoID << "\",\n";
                    outputFile << "\"ORF_start\": " << mem_lbit << ",\n";
                    outputFile << "\"ORF_end\": " << mem_rbit << ",\n";
                    //len不提供，让前端自己算
                    outputFile << "\"ORF_type\": \"" << (ORFtype == 2 ? "cORF" : "nBSJ-ORF") << "\",\n";
                    outputFile << "\"protein_name\":\"" << proteinName << "\",\n";
                    outputFile << "\"peptide\": \"" << line << "\"\n}";
                    if(lineCount == totalLineCount){
                        outputFile << "\n";
                    }
                    else{
                        outputFile << ",\n";
                    }
                }
                hasHit = false;
            }
        }
    }
    outputFile << "]";
    std::cout << "line count = " << lineCount << std::endl;
    std::cout << "total ORF = " << totalLineCount / 2 << std::endl;
    std::cout << "map hit = " << hit_ORF_count << std::endl;
    std::cout << "ORF overlong deleted: " << overlong_ORF_count << std::endl;
    std::cout << "cORF detected:" << cORF_count << std::endl;
    #pragma endregion
}
void get_human_IRES(){
    std::string filename("E:/workplace/Bio/circRNAdata/IRES_human.csv");
    // 打开文件
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << filename << std::endl;
        return;
    }
    std::string line;
    const int totalLineCount = 64574;
    int currentLineCount = 0;

    std::getline(file, line);//消去第一行表头行
    currentLineCount++;
    std::ofstream output_file("human_IRES_final.json");

    output_file << "[\n";
    //行内部遍历
    while (std::getline(file, line)){
        //汇报进度
        currentLineCount++;
        if(currentLineCount % 256 == 0){
            std::cout << "wait...completed" + std::to_string(currentLineCount * 100 / totalLineCount) << "%" << std::endl;
        }
        //这个文件的列项之间是用连续的空格分开的
        std::vector<std::string> lineStrVec(0);
        int l = 0, r = 0;

        for (; r < line.size();) //内部遍历line
        {
            //如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
            switch(line[r]){
                case '\"':
                    l = r;//l位于左侧引号
                    r++;
                    while(line[r]!='\"')r++;//r会停止在下引号位置
                    lineStrVec.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
                    r += 2;
                    l = r;
                    break;
                case ',':
                    lineStrVec.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }
        }
        //加上尾词
        lineStrVec.push_back(line.substr(l, r - l));
        output_file << "{\n";
        /*PK,isoform_ID,IRES_ID,score,identity,
         * alignment_length,mismatch,gapopen,evalue,bitscore,
         * iso_start,iso_end,IRES_start,IRES_end,IRES_seq_state,
         * iso_seq,IRES_seq,iso_seqlen,key,align_state*/


        if(lineStrVec[1][0] == '\"'){
            output_file << "\"isoform_ID\": " << lineStrVec[1] << ",\n";
        }
        else{
            output_file << "\"isoform_ID\": \"" << lineStrVec[1] << "\",\n";
        }
        output_file << "\"IRES_ID\": \"" << lineStrVec[2] << "\",\n";
        output_file << "\"score\": " << lineStrVec[3] << ",\n";
        output_file << "\"identity\": " << lineStrVec[4] << ",\n";
        output_file << "\"aln_len\": " << lineStrVec[5] << ",\n";
        output_file << "\"mismatch\": " << lineStrVec[6] << ",\n";
        output_file << "\"gapopen\": " << lineStrVec[7] << ",\n";
        output_file << "\"iso_start\": " << lineStrVec[10] << ",\n";
        output_file << "\"iso_end\": " << lineStrVec[11] << ",\n";
        output_file << "\"IRES_start\": " << lineStrVec[12] << ",\n";
        output_file << "\"IRES_end\": " << lineStrVec[13] << ",\n";
        output_file << "\"iso_seq\": \"" << lineStrVec[15] << "\",\n";
        output_file << "\"IRES_seq\": \"" << lineStrVec[16] << "\",\n";
        output_file << "\"align_state\": \"" << lineStrVec[19] << "\"\n";

        if(currentLineCount <= totalLineCount - 1){
            output_file << "},\n";
        }
        else{
            output_file << "}\n";
        }
    }
    output_file << "]";
    std::cout << "proc end, lineCount= " << currentLineCount << std::endl;
    // 关闭文件
    file.close();
    output_file.close();

}
void generate_mouse_bsjs_without_function(){
    std::string input_filename = inputDataDir + "isoform_mouse_alldata.txt";
    std::ifstream source_file(input_filename);
    std::string output_filename = "empty_function_mouse_bsj.json";
    std::ofstream  outputfile(output_filename);

    std::map<std::string, bsj_map_copy> bsj_isocount_map;
    long totalLineCount = 530400;
    int lineCount = 0;
    //  species: 'Mouse',
    //  BSJID: 'chr10|100069714|100076001|-', 来自vec[2]且需要map
    //  现在BSJ的格式还和人类不一样，我打算先不改了，貌似并不会造成什么严重后果
    //  genename: 'CPN1', 来自vec[9]
    //  geneid: 'ENSG00000120054',来自vec[8]
    //  isoformcount: 1, 来自map
    //  strand: '-', 来自vec[6]
    std::string line;
    std::getline(source_file, line);//跳过头行
    while(std::getline(source_file, line)){
        lineCount++;
        if(lineCount % 5304 == 1){
            std::cout << "generating map..." << lineCount*100 /totalLineCount << "%.\n";
        }
        std::vector<std::string> lineStrVec;
        split_line_into_strvec('\t', false, lineStrVec, line);
        //这里的ID用的是：应该换成|
        std::string transformed_BSJID = "";
        for(int i = 0; i < lineStrVec[2].size(); i++){
            if(lineStrVec[2][i] == ':'){
                lineStrVec[2][i] = '|';
                break;
            }
        }
        auto findResult = bsj_isocount_map.find(lineStrVec[2]);
        if(findResult != bsj_isocount_map.end()){
            findResult -> second.isoform_count ++;
            if(findResult -> second.strand == "")findResult -> second.strand = lineStrVec[6];
            if(findResult -> second.gene_id == "")findResult -> second.gene_id = lineStrVec[8];
            if(findResult -> second.gene_name == "")findResult -> second.gene_name = lineStrVec[9];
        }
        else{
            bsj_isocount_map.insert(std::make_pair(lineStrVec[2],
                                                   bsj_map_copy(lineStrVec[6],lineStrVec[8],lineStrVec[9],1)));
        }
    }
    lineCount = 0;
    totalLineCount = bsj_isocount_map.size();
    outputfile << "[\n";
    for(auto mapItem : bsj_isocount_map){
        lineCount++;
        if(lineCount% 4480 == 1){
            std::cout << "writing..." << lineCount * 100 / totalLineCount << "%\n";
        }
        outputfile << "{\n\"species\": \"Mouse\",\n";
        outputfile << "\"BSJID\": \"" << mapItem.first << "\",\n";
        outputfile << "\"genename\": \"" << mapItem.second.gene_name << "\",\n";
        outputfile << "\"geneid\": \"" << mapItem.second.gene_id << "\",\n";
        outputfile << "\"strand\": \"" << mapItem.second.strand << "\",\n";
        outputfile << "\"isoformcount\": " << mapItem.second.isoform_count << "\n}";
        if(lineCount == totalLineCount){
            outputfile << "\n";
        }
        else{
            outputfile << ",\n";
        }
    }
    outputfile << "]\n";
}
void reconstruct_RBP_report_from_all_batches(){
    std::map<int, std::string> coded_id2real_id_map;
    std::string map_file_name = outputTempDir + "map.out";
    std::ifstream map_file(map_file_name);
    std::string line;
    int lineCount = 0;
    int totalLineCount = 536940;
    while(std::getline(map_file, line)){
        lineCount ++;
        if(lineCount % 5369 == 1){
            std::cout << "building seq map... " << lineCount * 100 / totalLineCount << "%\n";
        }
        if(line.size() <= 1){
            continue;
        }
        //chr10:100044905-100045501Z100044905-100045501%1
        int breakPoint = 0;
        for(int i = 0; i < line.size(); i++){
            if(line[i] == 'Z'){
                line[i] = '|';
            }
            if(line[i] == '%'){
                breakPoint = i;
            }
        }
        std::string isoformID = substr_by_be(line, 0, breakPoint - 1);
        int seq_No = std::stoi(substr_by_be(line, breakPoint + 1, line.size() - 1));
        coded_id2real_id_map.insert(std::make_pair(seq_No, isoformID));
    }
    std::string outputFileName = outputTempDir + "RBP_output.out";
    std::ofstream outputFile(outputFileName);
    //int totalBatch = 157;
    int totalBatch = 157;
    outputFile << "[\n";
    for(int i = 1; i <= totalBatch; i++){
        std::stringstream ss("");
        ss << outputTempDir << "All_Predictions" << std::to_string(i) << ".txt";

        std::string report_file_name = ss.str();
        std::ifstream report_file(report_file_name);
        std::cout << "Processing filename = " << report_file_name << std::endl;
                  std::vector<std::string> report_article;
        std::string line;
        while(std::getline(report_file, line)){
            report_article.push_back(line);
        }
        std::cout << "processing batch" << i << ", article length = " << report_article.size() << " lines" << std::endl;
        reconstruct_RBP_report_from_batch(report_article, coded_id2real_id_map, outputFile);
        ss.str("");
    }
    outputFile << "]\n";
}
void correct_RBP_output(){
    std::string inputFileName = outputTempDir + "RBP_output_correct_v2.json";
    std::ifstream inputFile(inputFileName);

    std::string outputFileName = "RBP_output_correct_v3.json";
    std::ofstream outputFile(outputFileName);

    long long estimated_doc_length = 97106564;
    long long lineCount = 0;
    std::string line;
    while(std::getline(inputFile, line)){
        lineCount++;
        if(lineCount % 971065 == 1){
            std::cout << "correcting..." << lineCount * 100 / estimated_doc_length << "%." << std::endl;
        }
        if(line.size() >= 6 && substr_by_be(line, 0, 6) == "\"K-mers"){
            std::stringstream ss;
            int deletePos = 0;
            for(int i = 0; i < line.size(); i++){
                if(line[i] == '-'){
                    deletePos = i;
                    break;
                }
            }
            ss << substr_by_be(line, 0, deletePos - 1);
            ss << substr_by_be(line, deletePos + 1, line.size() - 1);
            outputFile << ss.str() << std::endl;
        }
        else{
            outputFile << line << std::endl;
        }
    }
}
void reconstruct_RBP_report_from_batch(std::vector<std::string> &article, std::map<int,std::string> &seqNoMap, std::ofstream &outputFile){
    std::string currentSeq;
    std::string currentProtein;
    for(int i = 0; i < article.size(); i++){
        if(article[i].size() <= 3)continue;
        if(substr_by_be(article[i], 0, 2) == "seq"){
            int seqNo = std::stoi(substr_by_be(article[i], 3, article[i].size() - 1));
            auto find_result = seqNoMap.find(seqNo);
            if(find_result != seqNoMap.end()){
                currentSeq = find_result->second;
            }
            else{
                std::cout << find_result->first << " != " << seqNo << std::endl;
            }
            continue;
        }
        if(substr_by_be(article[i], 0, 6) == "Report!"){
            //填写蛋白质
            int l = 0, r = 0;
            for(int k = 0; k < article[i].size(); k++){
                if(article[i][k] == ':'){
                    l = k + 2;
                    break;
                }
            }
            for(int k = l + 1; k < article[i].size(); k++){
                if(article[i][k] == '('){
                    r = k - 1;
                    break;
                }
            }
            currentProtein = substr_by_be(article[i], l, r);
            //倒查，找到上一个Sequence行
            std::vector<int> start_positions;
            std::vector<std::string> motifs;
            std::vector<std::string> kmers;
            //遍历每一个倒查行
            for(int j = i - 1; j >= 0 && substr_by_be(article[j], 0, 7) != "Sequence"; j--){
                std::string &lineRef = article[j];
                if(lineRef[0] >= '0' && lineRef[0] <= '9'){
                    int cutPos = 0;
                    for(int k = 0; k < lineRef.size(); k++){
                        if(lineRef[k] == ' '){
                            cutPos = k;
                            break;
                        }
                    }
                    start_positions.push_back(std::stoi(substr_by_be(lineRef, 0, cutPos - 1)));
                    for(int k = cutPos; k < lineRef.size(); k++){
                        if(lineRef[k] != ' '){
                            cutPos = k;
                            break;
                        }
                    }
                    int endPos = 0;
                    for(int k = cutPos; k < lineRef.size(); k++){
                        if(lineRef[k] == ' '){
                            endPos = k;
                            break;
                        }
                    }
                    motifs.push_back(substr_by_be(lineRef, cutPos + 1, endPos - 1));
                    for(int k = endPos; k < lineRef.size(); k++){
                        if(lineRef[k] != ' '){
                            cutPos = k;
                            break;
                        }
                    }
                    for(int k = cutPos; k < lineRef.size(); k++){
                        if(lineRef[k] == ' '){
                            endPos = k;
                            break;
                        }
                    }
                    kmers.push_back(substr_by_be(lineRef, cutPos + 1, endPos - 1));
                }
            }
            //写入JSON; 三个记录vector的size = 0的情况是不可能的，因为这里是Report检测行的if分支
            outputFile << "{\n";
            outputFile << "\"isoformID\": \"" << currentSeq << "\",\n";
            outputFile << "\"Protein:\": \"" << currentProtein <<"\",\n";
            outputFile << "\"BindingPos:\": [";
            for(int k = 0; k < start_positions.size() - 1; k++){
                outputFile << start_positions[k] << ", ";
            }
            outputFile << start_positions[start_positions.size() - 1] << "],\n";

            outputFile << "\"Motifs\": [";
            for(int k = 0; k < motifs.size() - 1; k++){
                outputFile << "\"" << motifs[k] << "\", ";
            }
            outputFile << "\"" << motifs[motifs.size() - 1] << "\"],\n";

            outputFile << "\"K-mers\": [";
            for(int k = 0; k < kmers.size() - 1; k++){
                outputFile << "\"" << kmers[k] << "\", ";
            }
            outputFile << "\"" << kmers[kmers.size() - 1] << "\"]\n},\n";
            /*
    * RBP分析程序有两项任务：1. 生成mouse_RBPs集合，这一集合的内容应该是这样的：
    * isoformID:xx
    * Protein: xx
    * Binding Positions: [12, 14, 15, 18]
    * Motifs: [uauau, uauuu, uauau, uaugu]
    * K-mers: [uauau, uauuu, uauau, uaugu]
    * 同一个isoformID将会有多个JSON，分别对应不同的Protein
    * */
        }
    }
}
void get_Mouse_IRES(){
    std::string input_filename = inputDataDir + "mmu.IRES.blat.out";
    std::ifstream source_file(input_filename);
    std::string output_filename = "mouseIRES.json";
    std::ofstream outputfile(output_filename);

    std::string line;
    outputfile << "[\n";
    int lineCount = 0;
    while(std::getline(source_file, line)){
        lineCount++;
        std::vector<std::string> lineStrVec;
        split_line_into_strvec('\t', false, lineStrVec, line);

        outputfile << "{\n";
        outputfile << "\"isoform_ID\": \"" << lineStrVec[1] << "\",\n";
        for(int i = 0; i < lineStrVec[0].size(); i++){
            if(lineStrVec[0][i] == '.'){
                lineStrVec[0] = substr_by_be(lineStrVec[0], 0, i - 1);
                break;
            }
        }
        outputfile << "\"IRES_ID\": \"" << lineStrVec[0] << "\",\n";
        outputfile << "\"identity\": " << lineStrVec[2] << ",\n";
        outputfile << "\"score\": " << lineStrVec[3] << ",\n";
        outputfile << "\"iso_start\": " << lineStrVec[6] << ",\n";
        outputfile << "\"iso_end\": " << lineStrVec[7] << ",\n";
        outputfile << "\"IRES_start\": " << lineStrVec[8] << ",\n";
        outputfile << "\"IRES_end\": " << lineStrVec[9] << "\n";
        outputfile << "}";
        if(lineCount == 1277){
            outputfile << "\n";
        }
        else{
            outputfile << ",\n";
        }
    }
    outputfile << "]";
}
void get_Mouse_miRNA(){
    std::string input_filename = inputDataDir + "isoform_mouse.fa.miranda.out.filt.cat.txt";
    std::ifstream source_file(input_filename);
    std::string output_filename = "mousemiRNA.json";
    std::ofstream outputfile(output_filename);

    std::string line;
    int lineCount = 0;
    int totalLineCount = 11881367;
    outputfile << "[\n";
    while(std::getline(source_file, line)){
        lineCount++;
        if(lineCount % 118813 == 1){
            std::cout << "waiting..." << lineCount * 100 / totalLineCount << "%" << std::endl;
        }
        std::vector<std::string> lineStrVec;
        split_line_into_strvec('\t', false, lineStrVec, line);
        outputfile << "{\n\"miRNAID\": \"" << lineStrVec[0] << "\",\n";
        outputfile << "\"isoformID\": \"" << lineStrVec[1] << "\", \n";
        outputfile << "\"score\":" << lineStrVec[2] << ",\n";
        outputfile << "\"energy\":" << lineStrVec[3] << ",\n";
        outputfile << "\"miRNAinteval\": \"" << lineStrVec[4] << "\",\n";
        outputfile << "\"isointeval\": \"" << lineStrVec[5] << "\",\n";
        outputfile << "\"alnlen\": " << lineStrVec[6] << ",\n";
        outputfile << "\"identity\": \"" << lineStrVec[7] << "\",\n";
        outputfile << "\"similarity\": \"" << lineStrVec[8] << "\"\n}";
        if(lineCount == totalLineCount){
            outputfile << "\n";
        }
        else{
            outputfile << ",\n";
        }
        //miRNAID           isoformID                                   score   energy  [miRNAinterval] isostart isoend identity similarity
        //mmu-let-7a-1-3p	chr10:12761224-12763101|12761224-12763101	179.00	-21.95	2          21	1117 1139	20	85.00%	90.00%
    }
    outputfile << "]\n";
    std::cout << lineCount << std::endl;

}
void get_Mouse_m6a_from_flat_text(){
    std::string input_filename = inputDataDir + "isoform_mouse.fa.miranda.out.filt.cat.txt";
    std::ifstream source_file(input_filename);
    std::string output_filename = "mousemiRNA.json";
    std::ofstream outputfile(output_filename);
}
void mainFun4MouseBSJ_Treatment(){
    //以下是对mouseBSJ.csv的处理程序
    //std::string filename = "/mnt/data/workdir/wangying/circRNA_nanopore/FLcirc_database/BSJ_mouse.csv";
    std::string filename = "your_file1.txt";
    // 打开文件
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件 " << filename << std::endl;
        return;
    }
    std::string line;
    long long lineCount = 0;
    long long totalLine = 115174;
    std::set<std::string> stored_Otherdb_str;

    while (std::getline(file, line)) {
        //进度提示器
        if (lineCount % 1150 == 1) {
            std::cout << "wait..." << lineCount / 115174 * 100 << "%" << std::endl;
        }
        //split by comma
        std::vector<std::string> splited(0);
        int l = 0;
        int r = 0;
        for (; r < line.size();) {
            //如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
            switch(line[r]){
                case '\"':
                    l = r;//l位于左侧引号
                    r++;
                    while(line[r]!='\"')r++;//r会停止在下引号位置
                    splited.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
                    r += 2;
                    l = r;
                    break;
                case ',':
                    splited.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }

        }
        splited.push_back(line.substr(l, r));//最后还会有一个忘却的record,此时l位于有效首位，r位于size-1

        stored_Otherdb_str.insert(splited[splited.size() - 9]);
    }
    std::cout << "proc end. " << std::endl;
    std::cout << "encountered Xref_DBs are:" << std::endl;
    for(auto str : stored_Otherdb_str){
        std::cout << str << std::endl;
    }
    // 关闭文件
    file.close();
}
//工具区
void long_isoID_encoder(){
    //原版数据必须经过这样几个步骤才能传递给BLAST用于建库：
    /*
     * 1. ID内|字符避讳
     * 2. ID超长，压缩为50个字符内
     * 这段程序完成这样一个过程：假定输入的fasta文件已经被ID——escape2charZ处理
     * 该段程序为所有长circRNAID建立编号映射，映射形成的map会写入encoded_isoID_map.txt
     * */
    std::string inputFastaFileName("./all_mouse_mm10_escape.fasta");
    std::ifstream inputFastaFile(inputFastaFileName);
    std::ofstream outputEncodedFasta("all_mouse_encoded.fasta");
    std::ofstream outputMap("map.out");

    if (!inputFastaFile.is_open()) {
        std::cerr << "无法打开文件 " << inputFastaFileName << std::endl;
        return;
    }
    std::string line;
    long long currentLineCount = 0;
    long long totalLineCount = 940350;

    std::string lastLine;
    while (std::getline(inputFastaFile, line)){
        //汇报进度
        currentLineCount++;
        if(currentLineCount % 9403 == 0){
            std::cout << "wait...completed" + std::to_string(currentLineCount * 100 / totalLineCount) << "%" << std::endl;
        }
        if(line.size() == 0){
            break;
        }
        if(line[0] == '>'){//ID抬头行处理
            std::string recorded_ID = line.substr(1, line.size()-1);//去掉第一个字符>就是需要的ID
            outputMap << recorded_ID << "%" << currentLineCount << std::endl;//这是map文件的一行，分隔符是%
            outputEncodedFasta << ">seq" << currentLineCount << std::endl;
        }
        else{
            //序列行处理
            outputEncodedFasta << line << std::endl;
        }
        lastLine = line;//更新lastline
        currentLineCount++;
    }

    // 关闭文件
    inputFastaFile.close();
    outputEncodedFasta.close();
    outputMap.close();

}
int ID_escape2charZ(){
    std::ifstream inputFile("./all_mouse.mm10.fa");

    if (!inputFile) {
        std::cerr << "无法打开文件。" << std::endl;
        return 1;
    }

    std::ofstream outputFile("all_mouse_mm10_escape.txt");

    if (!outputFile) {
        std::cerr << "无法创建输出文件。" << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        size_t pos = line.find('|');
        while (pos != std::string::npos) {
            line.replace(pos, 1, "Z");
            pos = line.find('|', pos + 1);
        }
        outputFile << line << std::endl;
    }
    std::cout << "替换完成。" << std::endl;
    inputFile.close();
    outputFile.close();
}
std::string substr_by_be(std::string& str, int b, int e){
    try{
        if(b > e || e >= str.size()){
            std::stringstream ss;
            ss << "substr_by_be_error, inputstr:" << str << " b = " << b << " e = " << e;
            throw std::runtime_error(ss.str());
        }
        return str.substr(b, e - b + 1);
    } catch (const std::runtime_error& e){
        std::cout << e.what() << std::endl;
    }
}
void split_line_into_strvec(char cutByChar, bool needQuote, std::vector<std::string> &lineStrVec, std::string &line){
    int l = 0, r = 0;
    for (; r < line.size();) //内部遍历line
    {
        if (needQuote && line[r] == '\"') {
            //如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
            l = r;//l位于左侧引号
            r++;
            while (line[r] != '\"')r++;//r会停止在下引号位置
            lineStrVec.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
            r += 2;
            l = r;
            continue;
        }
        if (line[r] == cutByChar) {
            lineStrVec.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
            r = r + 1;
            l = r;
            continue;
        }
        r++;
    }
    lineStrVec.push_back(line.substr(l, line.size() - l));
}
std::string interval_style_ID_convert2_array_style(std::string &ID, bool containBSJID, char strand){
    std::string result = interval_style_ID_convert2_array_style(ID, containBSJID);
    std::stringstream ss;
    ss << result;
    ss << "|" << strand;
    return ss.str();
}
std::string interval_style_ID_convert2_array_style(std::string &ID, bool containBSJID){
    std::stringstream ss;
    bool firstDivHasFound = false;
    int colonPos = 0;
    int divPos = 0;
    for(int i = 0; i < ID.size(); i++){
        if(ID[i] == ':'){
            ss << substr_by_be(ID, 0, i - 1);
            colonPos = i;
            continue;
        }
        if(ID[i] == '|'){
            firstDivHasFound = true;
            divPos = i;
            break;
        }
    }
    if(containBSJID){
        ss << substr_by_be(ID, colonPos, divPos);
    }
    else{
        ss << "|";
    }
    try{
        if(!firstDivHasFound){
            throw std::runtime_error("Input string Format Error:" + ID);
        }
        std::vector<std::string> seg_begins;
        std::vector<std::string> seg_ends;
        std::vector<int> cutPositions;
        for(int i = divPos + 1; i < ID.size(); i++){
            //intervalID例子：
            //chr10:100051846-100064146|100051846-100051959,100064084-100064146 //interval
            //chr10:100051846-100064146|100051846,100064084|100051959,100064146 //our
            //chr10|100051846,100064084|100051959,100064146 //our no BSJiD
            if(ID[i] == '-' || ID[i] == ','){
                cutPositions.push_back(i);
            }
        }
        cutPositions.push_back(ID.size());
        int l = divPos + 1;
        for(int k = 0; k < cutPositions.size()/2; k++){
            seg_begins.push_back(substr_by_be(ID, l, cutPositions[2 * k] - 1));
            l = cutPositions[2 * k] + 1;
            seg_ends.push_back(substr_by_be(ID, l, cutPositions[2 * k + 1] - 1));
            l = cutPositions[2 * k + 1] + 1;
        }
        for(int i = 0; i < seg_begins.size() - 1; i++){
            ss << seg_begins[i] << ",";
        }
        ss << seg_begins[seg_begins.size() - 1] << "|";

        for(int i = 0; i < seg_ends.size() - 1; i++){
            ss << seg_ends[i] << ",";
        }
        ss << seg_ends[seg_ends.size() - 1];
    }catch(const std::runtime_error& e){
        std::cout << e.what() << std::endl;
    }


    return ss.str();
}
int judgeORFvalid(int circLen, int l, int r){
    //当ORF长度超过一整个repeat长，这个ORF将被剔除；返回值为0
    //当ORF越过BSJ，这个ORF将被标识为cORF，返回值为2
    //否则，这是一个线性本和环形本共享的ORF，返回值为1
    int ll = l < r ? l : r;
    int rr = l < r ? r : l;
    if(rr - ll >= circLen)return 0;
    //总序列的下标应当是[0, 4 * circLen - 1],这当中的BSJ截断点在{circLen - 1; 2 * circLen - 1; 3 * circLen - 1; 4 * circLen - 1}
    if((ll + 1) / circLen != (rr + 1) / circLen)return 2;
    return 1;
}
void count_raw_RBP_batch_count(){
    long long totalCount = 0;
    for(int i = 1; i <= 116; i++){
        std::string RBPfilename = "RBPoutput_batch" + std::to_string(i) + ".out";
        std::ifstream inputFile(RBPfilename);
        std::string line;
        std::vector<std::string> lineStrVec;
        while(std::getline(inputFile, line)){
            split_line_into_strvec(' ', false, lineStrVec, line);
            totalCount += std::stoi(lineStrVec[2]);
        }
    }
    std::cout << totalCount << std::endl;
}
void automatic_scp_and_process_RBP_batch(){
    std::string scp_pass = "quaintSenator_1";
    std::system("cd E:/workplace/Bio/circRNAdata");
    for(int i = 1; i <= 23; i++){
        std::stringstream ss;
        ss  << "scp -i E:/workplace/Bio/circRNAdata/id_rsa " << "root@139.196.108.62:/usr/sf/resultbatch"
        << i << "/All_Predictions.txt ./All_Predictions" << i << ".txt";
        std::cout << "Transforming: resultbatch" << i << std::endl;
        int cmd_result = std::system(ss.str().c_str());
        if(cmd_result == 0){
            get_RBPdata_from_all_batch(i);
        }
    }
}

//绘图区
void get_Mouse_length_chromo_distribution_matlab_data(){
    std::string mouse_iso_bed_name = "E:/workplace/Bio/circRNAdata/isoform_mouse.bed";
    std::string outputFileName = "mouse_iso_distribution_matdat.dat";
    std::ifstream mouse_iso_file(mouse_iso_bed_name);
    std::ofstream outputFile(outputFileName);
    std::string line;
    int lineCount = 1;
    const long long totalLineCount = 536940;

    while(std::getline(mouse_iso_file, line)){
        if (lineCount % 5369 == 1) {
            std::cout << "reading file..." << lineCount * 100 / totalLineCount << "%\n";
        }
        lineCount++;
        std::vector<std::string> lineStrVec;

        int l = 0, r = 0;
        for (; r < line.size();) //内部遍历line
        {
            switch (line[r]) {
                case '\"'://如果发现了引号必须特殊处理，引号内部的逗号会扰乱split
                    l = r;//l位于左侧引号
                    r++;
                    while (line[r] != '\"')r++;//r会停止在下引号位置
                    lineStrVec.push_back(line.substr(l, r - l + 1)); // l to r (l, r - l + 1)
                    r += 2;
                    l = r;
                    break;
                case '\t':
                    lineStrVec.push_back(line.substr(l, r - l));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }
        }
        lineStrVec.push_back(line.substr(l, line.size() - l));
        std::string chromoNo = lineStrVec[0].substr(3, lineStrVec[0].size() - 3);
        //cut lineStrVec[10]
        std::string &lengthListRef = lineStrVec[10];
        l = r = 0;
        std::vector<int> lens;
        for (; r < lengthListRef.size();)
        {
            switch (lengthListRef[r]) {
                case ',':
                    lens.push_back(std::stoi(lengthListRef.substr(l, r - l)));// l to r - 1 (l, r- l)
                    r = r + 1;
                    l = r;
                    break;
                default:
                    r++;
                    break;
            }
        }
        lens.push_back(std::stoi(lengthListRef.substr(l, lengthListRef.size() - l)));
        int totalLen = 0;
        for(int i = 0; i < lens.size(); i++){
            totalLen += lens[i];
        }
        if(chromoNo[0] >= '0' && chromoNo[0] <= '9'){
            //数码染色体
        }
        else{
            if(chromoNo[0] == 'X') {
                chromoNo = "20";
            }
            else if(chromoNo[0] == 'Y'){
                chromoNo = "21";
            }
            else{
                chromoNo = "22";
            }
        }

        outputFile << chromoNo << "\t" << totalLen << std::endl;
    }

}
void get_Mouse_read_chromo_distribution_matlab_data(){
    //"TB","CC","HP","ST" 先readcount 后RPM
    /* x轴是四个分区，y轴是RPM画一张图，把y轴换成read count画一张图
     * x轴是四个分区，y轴是RPM，如果一个亚型TAU不为1，说明分别在不同区表达，那么这种点绘制成红色，否则绘制成蓝色
     * 还是做RPM的单一表
     */
    std::string mouse_exps_filename = outputTempDir + "allmouse_exp.bed";
    std::ifstream mouse_exps_file(mouse_exps_filename);
    std::string outputFileName = "spe_exp_matdat.dat";
    std::ofstream outputFile(outputFileName);

    std::string line;
    std::getline(mouse_exps_file, line);
    while(std::getline(mouse_exps_file, line)){
        std::vector<std::string> lineStrVec;
        split_line_into_strvec(' ', false, lineStrVec, line);
        std::vector<float> lineVecNum;
        for(int i = 1; i <= 9; i++){
            lineVecNum.push_back(std::stof(lineStrVec[i]));
        }
        //lineVecNum[0] = TAU
        for(int i = 5; i <= 8; i++){
            if(lineVecNum[i] != 0){
                //汇报RPM;一行有两个就都汇报
                int isTAU_1 = 0;
                isTAU_1 = (lineVecNum[0] == 1.0 ? 1 : 0);
                outputFile << (i - 4) << " " << lineVecNum[i] << " " << isTAU_1 << std::endl;
            }
        }
    }
}
void get_Mouse_chromo_RBP_chromosome_distribution_matlab_data(){
    std::string input_jsonFileName = outputTempDir + "RBP_output_correct_v3.json";
    std::ifstream file(input_jsonFileName);
    std::string outputDataFileName = outputTempDir + "RBP_heatmap_data4matlab.dat";
    std::ofstream outputData(outputDataFileName);
    std::map<std::string, int> protein2index;
    std::set<std::string> proteinNames;
    // 读取文件内容到字符串
    std::string line;
    bool insideOneJsonObj = false;
    long long totalLineCount = 97106562;
    long long lineCount = 0;
    std::vector<int> RBP_line(133, 0); //只使用[1 - 132]
    std::vector<std::vector<int>> matrixOfRBP_chr(23, RBP_line);//只使用 [1 - 22]
    int chrNo = 0;
    //这样使用matrixOfRBP_chr ---- matrixOfRBP_chr[chrNo][RBPNo]
    while(std::getline(file, line)){
        lineCount ++;
        if(lineCount % 971065 == 1){
            std::cout << "processing..." << lineCount * 100 / totalLineCount << "%" << std::endl;
        }
        /*if(substr_by_be(line, 0, 9) == "\"isoformID"){
        }*/
        if(line.size() <= 9)continue;
        if(substr_by_be(line, 0, 9) == "\"isoformID"){
            int start = 0, end = 0;
            for(int i = 9; i < line.size(); i++){
                if(line[i] == 'r'){
                    start = i + 1;
                    break;
                }
            }
            for(int i = start; i < line.size(); i++){
                if(line[i] == ':'){
                    end = i - 1;
                    break;
                }
            }
            if(line[start] >= '0' && line[start] <= '9'){
                chrNo = std::stoi(substr_by_be(line, start, end));
            }
            else{
                switch(line[start]){
                    case 'X':
                        chrNo = 20;
                        break;
                    case 'Y':
                        chrNo = 21;
                        break;
                    default:
                        chrNo = 22;
                        break;
                }
            }
        }
        else if(substr_by_be(line, 0, 7) == "\"Protein"){
            int begin = 0, end = 0;
            for(int i = 7; i < line.size(); i++){
                if(line[i] == ':'){
                    begin = i + 3;
                    break;
                }
            }
            for(int i = begin; i < line.size(); i++){
                if(line[i] == '\"'){
                    end = i - 1;
                    break;
                }
            }
            std::string id2insert = substr_by_be(line, begin, end);
            auto findResult = protein2index.find(id2insert);
            if(findResult == protein2index.end()){
                protein2index.insert(std::make_pair(id2insert, protein2index.size() + 1));
            }
            else{
                int RBPNo = findResult->second;
                matrixOfRBP_chr[chrNo][RBPNo]++;
            }
        }
    }
    for(auto pair : protein2index){
        std::cout << pair.first << "\t";
    }
    for(int i = 0; i < matrixOfRBP_chr.size(); i++){
        for(int j = 0; j < matrixOfRBP_chr[i].size(); j++){
            outputData << matrixOfRBP_chr[i][j] << " ";
        }
        outputData << std::endl;
    }
}
void report_our_identification_results(){
    //从获取的全亚型差异表达总表isoform_mouse_alldata.txt中提取每个亚型对每个组织差异表达矩阵
    //最终形成的差异表达矩阵的纵轴是所有亚型的ID，横轴是各个组织/细胞系，且每个组织/细胞系都有两个记录，分别对应isoc方法和CIRI方法
    std::string inputFileName = "E:/workplace/Bio/circRNAdata/isoform_mouse_alldata.txt";
    std::ifstream inputFile(inputFileName);

    // 打开文件
    if (!inputFile.is_open()) {
        std::cerr << "无法打开文件 " << inputFileName << std::endl;
        return;
    }
    const int totalLineCount = 536941;
    int lineCount = 1;
    std::string line;
    std::string headLine;

    //每个sample的总reads
    std::vector<float> totalReadVec(126, 0.0);
    std::vector<std::string> tableTabNames;
    //每个sample的总isoforms
    std::vector<int> countIsoVec(126, 0);
    std::getline(inputFile, headLine);
    int our_ciri_total_reads = 0;
    int our_ciri_total_isos = 0;
    int our_isocirc_total_reads = 0;
    int our_isocirc_total_isos = 0;
    int both_detected_isos = 0;
    split_line_into_strvec('\t', true, tableTabNames, headLine);

    while(std::getline(inputFile, line)){
        if(lineCount % 5360 == 1){
            std::cout << "counting..." << lineCount * 100 / totalLineCount << "%\n";
        }
        lineCount++;
        std::vector<std::string> lineStrVec;
        split_line_into_strvec('\t', false, lineStrVec, line);
        bool isoCircFound = false;
        bool cirilongFound = false;
        for(int i = 10; i <= 125; i++){
            float currentParam = std::stof(lineStrVec[i]);
            if(currentParam > 10000){
                std::cout << "BAD" << std::endl;
            }
            totalReadVec[i] += currentParam;//加和
            if(currentParam > 0.1){
                countIsoVec[i]++;
                if(i <= 101)cirilongFound = true;
                else isoCircFound = true;
            }
        }
        if(cirilongFound && isoCircFound)both_detected_isos++;
        else if(cirilongFound)our_ciri_total_isos++;
        else our_isocirc_total_isos++;
    }

    for(int i = 78; i <= 101; i++){//Ciri部分
        our_ciri_total_reads += countIsoVec[i];
    }
    for(int i = 102; i < countIsoVec.size(); i++){
        our_isocirc_total_reads += countIsoVec[i];
    }

    std::cout << "our_ciri_total_isos = " << our_ciri_total_isos << std::endl;
    std::cout << "our_isocirc_total_isos = " << our_isocirc_total_isos << std::endl;
    std::cout << "our_ciri_total_reads = " << our_ciri_total_reads << std::endl;
    std::cout << "our_isocirc_total_reads = " << our_isocirc_total_reads << std::endl;
    std::cout << "both: " << both_detected_isos << std::endl;
}
void get_Mouse_chromo_micro_chromosome_distribution_matlab_data(){
    std::string inputFileName = inputDataDir + "isoform_mouse.fa.miranda.out.filt.cat.txt";
    std::ifstream inputFile(inputFileName);
    std::string outputFileName = "matlab_micro_chromo_distribut.dat";
    std::ofstream outputFile(outputFileName);

    std::string line;
    std::string currentmiID = "mmu-let-7a-1-3p";
    std::string currentIsoID = "chr10:12761224-12763101|12761224-12763101";

    int totalLineCount = 12000000;
    int lineCount = 0;
    std::map<std::string, int> miID2index;
    std::vector<std::string> names;
    names.push_back(currentmiID);
    std::vector<int> init_line(23);//use 1-22
    std::vector<std::vector<int>> mi_chr_matrix(2000, init_line);
    int currentMi_Chr_count = 0;
    int currentChr = 10;
    while(std::getline(inputFile, line)){
        lineCount++;
        if(lineCount % 120000 == 1){
            std::cout << "estimating..." << lineCount * 100 / totalLineCount << "%\n";
        }
        std::vector<std::string> lineVec;
        split_line_into_strvec('\t', false, lineVec, line);
        if(lineVec[0] != currentmiID){
            //new miRNA
            miID2index.insert(std::make_pair(lineVec[0], miID2index.size() + 1));
            mi_chr_matrix[miID2index[currentmiID]][currentChr] = currentMi_Chr_count;
            currentmiID = lineVec[0];
        }
        else{
            int l = 0, r = 0;
            for(int i = 0; i < lineVec[1].size(); i++){
                if(lineVec[1][i] == 'r'){
                    l = i + 1;
                    break;
                }
            }
            for(int i = l; i < lineVec[1].size(); i++){
                if(lineVec[1][i] == ':'){
                    r = i - 1;
                    break;
                }
            }
            std::string tryGet_chrNo = substr_by_be(lineVec[1], l, r);
            int chrNo2Update = 0;
            if(tryGet_chrNo[0] >= '0' && tryGet_chrNo[0] <= '9'){
                chrNo2Update  = std::stoi(tryGet_chrNo);
            }
            else{
                switch(tryGet_chrNo[0]){
                    case 'X':
                        chrNo2Update = 20;
                        break;
                    case 'Y':
                        chrNo2Update = 21;
                        break;
                    default:
                        chrNo2Update = 22;
                        break;
                }
            }
            if(chrNo2Update != currentChr){
                //new chr
                mi_chr_matrix[miID2index[currentmiID]][currentChr] = currentMi_Chr_count;
                if(currentMi_Chr_count > 4000){
                    std::cout << currentmiID << "! chr"<< currentChr << ";" << currentMi_Chr_count << std::endl;
                }
                currentMi_Chr_count = 1;
                currentChr = chrNo2Update;

            }
            else{
                //same chr_mi
                currentMi_Chr_count++;
            }
        }
    }
    mi_chr_matrix[miID2index[currentmiID]][currentChr] = currentMi_Chr_count;
    /*for(int i = 1; i <= 1937; i++){
        for(int j = 1; j < mi_chr_matrix[i].size(); j++){
            outputFile << mi_chr_matrix[i][j] << " ";
            for(int k = 1; k <= 83; k++){
                outputFile << 0 << " ";
            }
        }
        outputFile << std::endl;
    }*/
    for(int i = 1; i <= 1937; i++){
        for(int j = 1; j < mi_chr_matrix[i].size(); j++){
            outputFile << mi_chr_matrix[i][j] << " ";
            for(int k = 1; k <= 83; k++){
                outputFile << 0 << " ";
            }
        }
        outputFile << std::endl;
    }
}