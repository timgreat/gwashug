import java.io.*;
import java.util.*;
public class Function {
    String dataFile;
    String workDir;

    String ancestry;
    String trait;
    String mrpvalue;
    String samestudy;
    String type;
    String usemodel1;
    String usemodel2;
    String pph4;
    String cases;

    HashMap<String, String> params;

    String preprocess;
    public Function(){
        dataFile = "Data_file";
        workDir = "Work_dir";

        ancestry = "Ancestry";
        trait="Trait";
        mrpvalue="MRPvalue";
        samestudy="SameStudy";
        type="Trait_type";
        usemodel1="Tissue1";
        usemodel2="Tissue2";
        pph4="PPH4";
        cases="Number_of_Cases";

        String[] keys = new String[]{dataFile, workDir, ancestry,trait, mrpvalue,samestudy,type,usemodel1,usemodel2,pph4,cases};
        params = new HashMap<>();
        //init
        for (String k : keys)
            params.put(k, "");

        preprocess = "/root/anaconda3/bin/python /home/chencao/preprocess_gwas/preprocess.py ";
    }
    public void loadParams(String p){

        //init
        params.put(ancestry,"EUR");
        params.put(trait,"trait");
        params.put(mrpvalue,"1E-5");
        params.put(samestudy,"True");
        params.put(type,"Continuous");
        params.put(usemodel1,"GTEx_Whole_Blood");
        params.put(usemodel2,"GTEx_Whole_Blood");
        params.put(pph4,"0.5");
        params.put(cases,"non");

        BufferedReader bf = null;
        try {
            if (!new File(p).exists()){
                System.out.println(p + " doesn't exist!");
                System.exit(0);
            }

            bf = new BufferedReader(new FileReader(p));
            String line;
            while ((line = bf.readLine()) != null) {
                String[] pairs = line.split("\\s+");
                if (params.containsKey(pairs[0])){
                    if(pairs.length>=2){
                        String value=pairs[1];
                        for(int i=2;i<pairs.length;i++)
                            value=value+"_"+pairs[i];
                        params.put(pairs[0], value);
                    }

                }
            }
        }catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(p + " loaded.");
    }
    public String checkMRPvalue(String summ,String mrpvalue) throws IOException {
        double pva=Double.parseDouble(mrpvalue);
        BufferedReader bf=new BufferedReader(new FileReader(summ));
        String line=bf.readLine();
        double minPva=1.0;
        while ((line=bf.readLine())!=null){
            String[] snp=line.split("\\s+");
            double now=Double.parseDouble(snp[10]);
            if(minPva>now){

                minPva=now;
            }
        }
        if(pva>minPva)
            return mrpvalue;
        else
            return String.valueOf(minPva);
    }

//    public void pngList(String path) throws IOException {
//        File dir=new File(path);
//        if(dir.exists()){
//            File[] files=dir.listFiles();
//            FileWriter writer = new FileWriter(path+"/pictures.txt");
//            for(File file:files){
//                String fileName=file.getName();
//                if(fileName.contains(".png"))
//                    writer.write(fileName+"\n");
//            }
//            writer.close();
//        }
//    }
//    public void txtList(String path) throws IOException {
//        File dir=new File(path);
//        if(dir.exists()){
//            File[] files=dir.listFiles();
//            FileWriter writer = new FileWriter(path+"/download.txt");
//            for(File file:files){
//                String fileName=file.getName();
//                if(fileName.contains(".txt"))
//                    writer.write(fileName+"\n");
//            }
//            writer.close();
//        }
//    }
//    public void methodList(String path,String tissue1,String tissue2) throws IOException {
//        File dir=new File(path);
//        if(dir.exists()){
//            File[] files=dir.listFiles();
//            FileWriter writer = new FileWriter(path+"/method.txt");
//            writer.write(tissue1+","+tissue2+"\n");
//            String gene="";
//            int count=0;
//            for(File file:files){
//                String fileName=file.getName();
//                if(fileName.contains("summ1.png")){
//                    String[] array=fileName.split("_");
//                    if(count==0){
//                        gene=gene+array[array.length-2];
//                    }else {
//                        gene=gene+","+array[array.length-2];
//                    }
//                    count++;
//                }
//            }
//            writer.write(gene+"\n");
//            writer.write("summ1,summ2,eqtl1,eqtl2\n");
//            writer.close();
//        }
//    }
    public void methodList(String path,String tissue1,String tissue2) throws IOException {
        File dir=new File(path);
        if(dir.exists()){
            File[] files=dir.listFiles();
            FileWriter writer = new FileWriter(path+"/method.txt");
            writer.write(tissue1+","+tissue2+"\n");
            String gene1="",gene2="";
            int count1=0,count2=0;
            for(File file:files){
                String fileName=file.getName();
                if(fileName.contains(tissue1+"_summ1.png")){
                    String[] array1=fileName.split("_");
                    if(count1==0){
                        gene1=gene1+array1[array1.length-2];
                    }else {
                        gene1=gene1+","+array1[array1.length-2];
                    }
                    count1++;
                } else if(fileName.contains(tissue2+"_summ1.png")){
                    String[] array2=fileName.split("_");
                    if(count2==0){
                        gene2=gene2+array2[array2.length-2];
                    }else {
                        gene2=gene2+","+array2[array2.length-2];
                    }
                    count2++;
                }
            }
            writer.write(gene1+"\n");
            writer.write(gene2+"\n");
            writer.write("summ1,summ2,eqtl1,eqtl2\n");
            writer.close();
        }
    }

    //trait level
    //MR
    public void MRExec(){
        String summary1 = params.get(workDir) + "/trait1.txt";
        String mr_dir = params.get(workDir) + "/MR";
        String summary2 = params.get(workDir) + "/trait2.txt";
        try {
            if (!new File(summary1).exists() || !new File(summary2).exists()){
                new ProcessBuilder("sh", "-c", preprocess + params.get(workDir)).start().waitFor();
                System.out.println("preprocess finished!!!");
            }
            File dir = new File(mr_dir);
            if (!dir.exists())
                dir.mkdir();

            String summary_info_path = mr_dir + "/summary_info.txt";
            FileWriter writer = new FileWriter(summary_info_path);
            writer.write("GWAS1:\t" + summary1 + "\n");
            writer.write("Pop1:\t"+params.get(ancestry)+"\n");
            writer.write("Trait1:\t"+params.get(trait).replace(" ","_")+"\n");

            loadParams(params.get(workDir) + "/MR.xparam");
            writer.write("GWAS2:\t" + summary2 + "\n");
            writer.write("Pop2:\t"+params.get(ancestry)+"\n");
            writer.write("Trait2:\t"+params.get(trait).replace(" ","_")+"\n");

            String mrp=checkMRPvalue(summary1,params.get(mrpvalue));
            writer.write("MRPvalue:\t"+mrp+"\n");
            writer.close();

            String sh = "/root/anaconda3/envs/gwashug/bin/Rscript /root/cogwas/trait_level/MR.R --info " + summary_info_path;
            System.out.println(sh);

            ProcessBuilder processBuilder = new ProcessBuilder("sh", "-c", sh);
            Process process = processBuilder.start();
            InputStream inputStream = process.getInputStream();
            Scanner scanner = new Scanner(inputStream);
            while (scanner.hasNextLine()) {
                System.out.println(scanner.nextLine());
            }
            process.waitFor();

            System.out.println("MR FINISHED.");
        } catch (Exception e) {
            System.out.println("MR ERROR!!!");
            e.printStackTrace();
        }
    }
    //GECKO
    public void GECKOExec(){
        String summary1 = params.get(workDir) + "/trait1.txt";
        String prs_dir = params.get(workDir) + "/GECKO";
        String summary2 = params.get(workDir) + "/trait2.txt";
        try {
            if (!new File(summary1).exists() || !new File(summary2).exists()){
                new ProcessBuilder("sh", "-c", preprocess + params.get(workDir)).start().waitFor();
                System.out.println("preprocess finished!!!");
            }
            File dir = new File(prs_dir);
            if (!dir.exists())
                dir.mkdir();

            String summary_info_path = prs_dir + "/summary_info.txt";
            FileWriter writer = new FileWriter(summary_info_path);
            writer.write("GWAS1:\t" + summary1 + "\n");
            writer.write("Pop1:\t"+params.get(ancestry)+"\n");
            //writer.write("Trait1:\t"+params.get(trait)+"\n");

            loadParams(params.get(workDir) + "/GECKO.xparam");
            writer.write("GWAS2:\t" + summary2 + "\n");
            writer.write("Pop2:\t"+params.get(ancestry)+"\n");
            //writer.write("Trait2:\t"+params.get(trait)+"\n");

            if(params.get(samestudy).toLowerCase().equals("false"))
                writer.write("SameStudy:\tF\n");
            else
                writer.write("SameStudy:\tT\n");

            writer.close();

            String sh = "/root/anaconda3/envs/gwashug/bin/Rscript /root/cogwas/trait_level/GECKO.R  --info " + summary_info_path;
            System.out.println(sh);

            ProcessBuilder processBuilder = new ProcessBuilder("sh", "-c", sh);
            Process process = processBuilder.start();
            InputStream inputStream = process.getInputStream();
            Scanner scanner = new Scanner(inputStream);
            while (scanner.hasNextLine()) {
                System.out.println(scanner.nextLine());
            }
            process.waitFor();

            System.out.println("GECKO FINISHED.");
        } catch (Exception e) {
            System.out.println("GECKO ERROR!!!");
            e.printStackTrace();
        }
    }
    //LDSC
    public void LDSCExec(){
        String summary1 = params.get(workDir) + "/trait1.txt";
        String prs_dir = params.get(workDir) + "/LDSC";
        String summary2 = params.get(workDir) + "/trait2.txt";
        try {
            if (!new File(summary1).exists() || !new File(summary2).exists()){
                new ProcessBuilder("sh", "-c", preprocess + params.get(workDir)).start().waitFor();
                System.out.println("preprocess finished!!!");
            }
            File dir = new File(prs_dir);
            if (!dir.exists())
                dir.mkdir();

            String summary_info_path = prs_dir + "/summary_info.txt";
//            FileWriter writer = new FileWriter(summary_info_path);
//            writer.write("GWAS1:\t" + summary1 + "\n");
//            writer.write("Pop1:\t"+params.get(ancestry)+"\n");
            //writer.write("Trait1:\t"+params.get(trait)+"\n");
            String pop1=params.get(ancestry);

            loadParams(params.get(workDir) + "/LDSC.xparam");
//            writer.write("GWAS2:\t" + summary2 + "\n");
//            writer.write("Pop2:\t"+params.get(ancestry)+"\n");
            //writer.write("Trait2:\t"+params.get(trait)+"\n");
            String pop2=params.get(ancestry);

//            if(params.get(samestudy).toLowerCase().equals("false"))
//                writer.write("SameStudy:\tF\n");
//            else
//                writer.write("SameStudy:\tT\n");
//
//            writer.close();

            String sh = "/root/anaconda3/envs/gwashug/bin/Rscript /root/cogwas/trait_level/ldsc.R  --info " + summary_info_path;
            System.out.println(sh);

            ProcessBuilder processBuilder = new ProcessBuilder("sh", "-c", sh);
            Process process = processBuilder.start();
            InputStream inputStream = process.getInputStream();
            Scanner scanner = new Scanner(inputStream);
            while (scanner.hasNextLine()) {
                System.out.println(scanner.nextLine());
            }
            process.waitFor();

            System.out.println("LDSC FINISHED.");
        } catch (Exception e) {
            System.out.println("LDSC ERROR!!!");
            e.printStackTrace();
        }
    }
    //LAVA
    public void LAVAExec(){
        String summary1 = params.get(workDir) + "/trait1.txt";
        String prs_dir = params.get(workDir) + "/LAVA";
        String summary2 = params.get(workDir) + "/trait2.txt";
        try {
            if (!new File(summary1).exists() || !new File(summary2).exists()){
                new ProcessBuilder("sh", "-c", preprocess + params.get(workDir)).start().waitFor();
            }
            File dir = new File(prs_dir);
            if (!dir.exists())
                dir.mkdir();

            String summary_info_path = prs_dir + "/summary_info.txt";
            FileWriter writer = new FileWriter(summary_info_path);
            writer.write("GWAS1:\t" + summary1 + "\n");
            writer.write("Pop1:\t"+params.get(ancestry)+"\n");
            writer.write("Trait1:\t"+params.get(trait).replace(" ","_")+"\n");
            if (params.get(type).toLowerCase().equals("binary"))
                writer.write("Type1:\tcc\n");
            else
                writer.write("Type1:\tquant\n");
            String case1=params.get(cases).replace(",","");

            loadParams(params.get(workDir) + "/LAVA.xparam");
            writer.write("GWAS2:\t" + summary2 + "\n");
            writer.write("Pop2:\t"+params.get(ancestry)+"\n");
            writer.write("Trait2:\t"+params.get(trait).replace(" ","_")+"\n");
            if (params.get(type).toLowerCase().equals("binary"))
                writer.write("Type2:\tcc\n");
            else
                writer.write("Type2:\tquant\n");

            String case2=params.get(cases).replace(",","");
            if(!case1.equals("non")){
                writer.write("Case1:\t"+case1+"\n");
            }
            if(!case2.equals("non")){
                writer.write("Case2:\t"+case2+"\n");
            }
            writer.close();

            String sh = "/root/anaconda3/envs/gwashug/bin/Rscript /root/cogwas/trait_level/LAVA.R --info " + summary_info_path;
            System.out.println(sh);

            ProcessBuilder processBuilder = new ProcessBuilder("sh", "-c", sh);
            Process process = processBuilder.start();
            InputStream inputStream = process.getInputStream();
            Scanner scanner = new Scanner(inputStream);
            while (scanner.hasNextLine()) {
                System.out.println(scanner.nextLine());
            }
            process.waitFor();

            System.out.println("LAVA FINISHED.");
        } catch (Exception e) {
            System.out.println("LAVA ERROR!!!");
            e.printStackTrace();
        }
    }
    //COLOC
    public void COLOCExec(){
        String summary1 = params.get(workDir) + "/trait1.txt";
        String prs_dir = params.get(workDir) + "/COLOC";
        String summary2 = params.get(workDir) + "/trait2.txt";
        try {
            if (!new File(summary1).exists() || !new File(summary2).exists()){
                new ProcessBuilder("sh", "-c", preprocess + params.get(workDir)).start().waitFor();
            }
            File dir = new File(prs_dir);
            if (!dir.exists())
                dir.mkdir();

            String summary_info_path = prs_dir + "/summary_info.txt";
            FileWriter writer = new FileWriter(summary_info_path);
            writer.write("GWAS1:\t" + summary1 + "\n");
            writer.write("Pop1:\t"+params.get(ancestry)+"\n");
            writer.write("Trait1:\t"+params.get(trait).replace(" ","_")+"\n");
            if (params.get(type).toLowerCase().equals("binary"))
                writer.write("Type1:\tcc\n");
            else
                writer.write("Type1:\tquant\n");

            loadParams(params.get(workDir) + "/COLOC.xparam");
            writer.write("GWAS2:\t" + summary2 + "\n");
            writer.write("Pop2:\t"+params.get(ancestry)+"\n");
            writer.write("Trait2:\t"+params.get(trait).replace(" ","_")+"\n");
            if(params.get(type).toLowerCase().equals("binary"))
                writer.write("Type2:\tcc\n");
            else
                writer.write("Type2:\tquant\n");
            writer.write("PPH4:\t"+params.get(pph4)+"\n");

            writer.close();

            String sh = "/root/anaconda3/envs/gwashug/bin/Rscript /root/cogwas/trait_level/coloc.R --info " + summary_info_path;
            System.out.println(sh);

            ProcessBuilder processBuilder = new ProcessBuilder("sh", "-c", sh);
            Process process = processBuilder.start();
            InputStream inputStream = process.getInputStream();
            Scanner scanner = new Scanner(inputStream);
            while (scanner.hasNextLine()) {
                System.out.println(scanner.nextLine());
            }
            process.waitFor();
            System.out.println("COLOC FINISHED.");
        } catch (Exception e) {
            System.out.println("COLOC ERROR!!!");
            e.printStackTrace();
        }
    }

    //molecular level
    //COLOC_EQTL
    public void COLOCEQTLExec(){
        String summary1 = params.get(workDir) + "/trait1.txt";
        String prs_dir = params.get(workDir) + "/COLOCEQTL";
        String summary2 = params.get(workDir) + "/trait2.txt";
        try {
            if (!new File(summary1).exists() || !new File(summary2).exists()){
                new ProcessBuilder("sh", "-c", preprocess + params.get(workDir)).start().waitFor();
            }
            File dir = new File(prs_dir);
            if (!dir.exists())
                dir.mkdir();

            String summary_info_path = prs_dir + "/summary_info.txt";
            FileWriter writer = new FileWriter(summary_info_path);
            writer.write("GWAS1:\t" + summary1 + "\n");
            if(params.get(type).toLowerCase().equals("binary"))
                writer.write("Type1:\tcc\n");
            else
                writer.write("Type1:\tquant\n");

            loadParams(params.get(workDir) + "/COLOCEQTL.xparam");
            writer.write("GWAS2:\t" + summary2 + "\n");
            if(params.get(type).toLowerCase().equals("binary"))
                writer.write("Type2:\tcc\n");
            else
                writer.write("Type2:\tquant\n");
            writer.write("Use_model1:\t"+params.get(usemodel1)+"\n");
            writer.write("Use_model2:\t"+params.get(usemodel2)+"\n");
            writer.write("PPH4:\t"+params.get(pph4)+"\n");
            writer.close();

            String sh = "/root/anaconda3/envs/gwashug/bin/Rscript /root/cogwas/molecular_level/coloc_eqtl.R --info " + summary_info_path;
            System.out.println(sh);

            ProcessBuilder processBuilder = new ProcessBuilder("sh", "-c", sh);
            Process process = processBuilder.start();
            InputStream inputStream = process.getInputStream();
            Scanner scanner = new Scanner(inputStream);
            while (scanner.hasNextLine()) {
                System.out.println(scanner.nextLine());
            }
            process.waitFor();

            methodList(prs_dir+"/output",params.get(usemodel1),params.get(usemodel2));

            System.out.println("COLOC_EQTL FINISHED.");
        } catch (Exception e) {
            System.out.println("COLOC_EQTL ERROR!!!");
            e.printStackTrace();
        }
    }
    //SMR
    public void SMRExec(){
        String summary1 = params.get(workDir) + "/trait1.txt";
        String prs_dir = params.get(workDir) + "/SMR";
        String summary2 = params.get(workDir) + "/trait2.txt";
        try {
            if (!new File(summary1).exists() || !new File(summary2).exists()){
                new ProcessBuilder("sh", "-c", preprocess + params.get(workDir)).start().waitFor();
            }
            File dir = new File(prs_dir);
            if (!dir.exists())
                dir.mkdir();

            String summary_info_path = prs_dir + "/summary_info.txt";
            FileWriter writer = new FileWriter(summary_info_path);
            writer.write("GWAS1:\t" + summary1 + "\n");

            loadParams(params.get(workDir) + "/SMR.xparam");
            writer.write("GWAS2:\t" + summary2 + "\n");
            writer.write("Use_model1:\t"+params.get(usemodel1)+"\n");
            writer.write("Use_model2:\t"+params.get(usemodel2)+"\n");
            writer.close();

            String sh = "/root/anaconda3/envs/gwashug/bin/Rscript /root/cogwas/molecular_level/SMR.R --info " + summary_info_path;
            System.out.println(sh);

            ProcessBuilder processBuilder = new ProcessBuilder("sh", "-c", sh);
            Process process = processBuilder.start();
            InputStream inputStream = process.getInputStream();
            Scanner scanner = new Scanner(inputStream);
            while (scanner.hasNextLine()) {
                System.out.println(scanner.nextLine());
            }
            process.waitFor();

            FileWriter smr_writer = new FileWriter(prs_dir+"/output/method.txt");
            smr_writer.write(params.get(usemodel1)+","+params.get(usemodel2)+"\n");
            smr_writer.close();

            System.out.println("SMR FINISHED.");
        } catch (Exception e) {
            System.out.println("SMR ERROR!!!");
            e.printStackTrace();
        }
    }
    public void INTACTExec(){
        String summary1 = params.get(workDir) + "/trait1.txt";
        String prs_dir = params.get(workDir) + "/INTACT";
        String summary2 = params.get(workDir) + "/trait2.txt";
        try {
            if (!new File(summary1).exists() || !new File(summary2).exists()){
                new ProcessBuilder("sh", "-c", preprocess + params.get(workDir)).start().waitFor();
            }
            File dir = new File(prs_dir);
            if (!dir.exists())
                dir.mkdir();

            String summary_info_path = prs_dir + "/summary_info.txt";
            FileWriter writer = new FileWriter(summary_info_path);
            writer.write("GWAS1:\t" + summary1 + "\n");

            loadParams(params.get(workDir) + "/INTACT.xparam");
            writer.write("GWAS2:\t" + summary2 + "\n");
            writer.write("Use_model1:\t"+params.get(usemodel1)+"\n");
            writer.write("Use_model2:\t"+params.get(usemodel2)+"\n");
            writer.close();

            String sh = "/root/anaconda3/envs/gwashug/bin/Rscript /root/cogwas/molecular_level/INTACT.R --info " + summary_info_path;
            System.out.println(sh);

            ProcessBuilder processBuilder = new ProcessBuilder("sh", "-c", sh);
            Process process = processBuilder.start();
            InputStream inputStream = process.getInputStream();
            Scanner scanner = new Scanner(inputStream);
            while (scanner.hasNextLine()) {
                System.out.println(scanner.nextLine());
            }
            process.waitFor();

            FileWriter intact_writer = new FileWriter(prs_dir+"/output/method.txt");
            intact_writer.write(params.get(usemodel1)+","+params.get(usemodel2)+"\n");
            intact_writer.close();

            System.out.println("INTACT FINISHED.");
        } catch (Exception e) {
            System.out.println("INTACT ERROR!!!");
            e.printStackTrace();
        }
    }
    public void fbp(String path){
        long s, e;
        //MR
        if (new File(path + "/MR.param").exists()
                && new File(path + "/MR.xparam").exists()){
            s = System.currentTimeMillis();
            loadParams(path + "/MR.param");
            MRExec();
            e = System.currentTimeMillis();
            System.out.println("time consuming:" + (e - s) / 60000 + "min");
            System.out.println("-------------------------------");
        }
        //GECKO
        if (new File(path + "/GECKO.param").exists()
                && new File(path + "/GECKO.xparam").exists()){
            s = System.currentTimeMillis();
            loadParams(path + "/GECKO.param");
            GECKOExec();
            e = System.currentTimeMillis();
            System.out.println("time consuming:" + (e - s) / 60000 + "min");
            System.out.println("-------------------------------");
        }
        //COLOC
        if (new File(path + "/COLOC.param").exists()
                && new File(path + "/COLOC.xparam").exists()){
            s = System.currentTimeMillis();
            loadParams(path + "/COLOC.param");
            COLOCExec();
            e = System.currentTimeMillis();
            System.out.println("time consuming:" + (e - s) / 60000 + "min");
            System.out.println("-------------------------------");
        }
        //LAVA
        if (new File(path + "/LAVA.param").exists()
                && new File(path + "/LAVA.xparam").exists()){
            s = System.currentTimeMillis();
            loadParams(path + "/LAVA.param");
            LAVAExec();
            e = System.currentTimeMillis();
            System.out.println("time consuming:" + (e - s) / 60000 + "min");
            System.out.println("-------------------------------");
        }

        //COLOC_EQTL
        if (new File(path + "/COLOCEQTL.param").exists()
                && new File(path + "/COLOCEQTL.xparam").exists()){
            s = System.currentTimeMillis();
            loadParams(path + "/COLOCEQTL.param");
            COLOCEQTLExec();
            e = System.currentTimeMillis();
            System.out.println("time consuming:" + (e - s) / 60000 + "min");
            System.out.println("-------------------------------");
        }
        //SMR
        if (new File(path + "/SMR.param").exists()
                && new File(path + "/SMR.xparam").exists()){
            s = System.currentTimeMillis();
            loadParams(path + "/SMR.param");
            SMRExec();
            e = System.currentTimeMillis();
            System.out.println("time consuming:" + (e - s) / 60000 + "min");
            System.out.println("-------------------------------");
        }
        //INTACT
        if (new File(path + "/INTACT.param").exists()
                && new File(path + "/INTACT.xparam").exists()){
            s = System.currentTimeMillis();
            loadParams(path + "/INTACT.param");
            INTACTExec();
            e = System.currentTimeMillis();
            System.out.println("time consuming:" + (e - s) / 60000 + "min");
            System.out.println("-------------------------------");
        }
    }

}
