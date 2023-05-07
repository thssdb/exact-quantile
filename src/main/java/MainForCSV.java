import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.Random;

public class MainForCSV {
    static final int N = 100000000;

    public static void main(String[] args) throws IOException {
        Spacecraft_Thruster_Firing_Tests_Dataset();
    }

    public static void splice() throws IOException {
        int txt_id = 2;
        String filename = "Admiral_Schmidt-ems1-group1.csv";
        String folder = "D:\\Study\\Lab\\iotdb\\real-world data\\中船\\合并后数据";
        String filepath = folder + "\\" + filename;
        File file = new File(filepath);
        BufferedInputStream fis = null;
        fis = new BufferedInputStream(new FileInputStream(file));
        BufferedReader reader = new BufferedReader(new InputStreamReader(fis, StandardCharsets.UTF_8), 5 * 1024 * 1024);// 用5M的缓冲读取文本文件

        File oFile = new File(folder + "\\tmp_" + txt_id + ".txt");
        BufferedOutputStream fos = new BufferedOutputStream(new FileOutputStream(oFile));
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(fos, StandardCharsets.UTF_8), 5 * 1024 * 1024);
        writer.write("from: " + folder + " " + filename + "\n");

        String line = "";
        int id = 0;
        long tmp = 0;
        double[] lastVal = new double[23333];
        line = reader.readLine();//ignore first line
        while ((line = reader.readLine()) != null) {
            //TODO: write your business
            id++;
//            System.out.println("\t\t"+line);
            String[] strArr = line.split(",");
            double tmpV;
            for (int i = 1; i < strArr.length; i++)
                if (strArr[i].length() > 0) {
                    try {
                        tmpV = Double.parseDouble(strArr[i]);
                    } catch (NumberFormatException e) {
                        continue;
                    }
//                System.out.println("\t\t"+tmpV);
                    if (tmpV != lastVal[i]) {
                        lastVal[i] = tmpV;
                        tmp++;

                        writer.write(tmpV + "\n");
                        if (tmp >= N) break;
//                    if(tmp<=50)
//                    System.out.println("\t\t"+tmpV);
                    }
                }
            if (tmp >= N) break;
//            if(id>=50)break;
        }
        reader.close();
        fis.close();
        System.out.println("?\t\t" + tmp);
        writer.close();
        fos.close();
    }

    public static void ExtractInFolderN() throws IOException {
        final int SAME_LIMIT = 30;
        String folder = "D:\\Study\\Lab\\iotdb\\real-world data\\Kaggle\\Household appliances power consumption";
        File folderFile = new File(folder);
        File[] files = folderFile.listFiles();
        assert files != null;
        int txt_id = 0;

        if (!(new File(folder + "\\N" + N)).exists()) (new File(folder + "\\N" + N)).mkdir();

        File oFile = new File(folder + "\\N" + N + "\\tmp_" + txt_id + ".txt");
        BufferedOutputStream fos = new BufferedOutputStream(new FileOutputStream(oFile));
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(fos, StandardCharsets.UTF_8), 5 * 1024 * 1024);
        writer.write("from: " + folder + "\n");
        long tmp = 0;
        int same = 0;

        for (File file : files)
            if (!file.isDirectory()) {
                System.out.println("\t\t\t" + file.getName());
                if (!file.getName().contains(".csv")) continue;
                BufferedInputStream fis = null;
                fis = new BufferedInputStream(new FileInputStream(file));
                BufferedReader reader = new BufferedReader(new InputStreamReader(fis, StandardCharsets.UTF_8), 5 * 1024 * 1024);// 用5M的缓冲读取文本文件

                String line = "";
                int id = 0;
                double lastVal = -1;
                line = reader.readLine();//ignore first line
                while ((line = reader.readLine()) != null) {
                    //TODO: write your business
                    id++;
//            System.out.println("\t\t"+line);
                    String[] strArr = line.split(",");
                    double tmpV;
                    for (int i = 1; i < strArr.length; i++)
                        if (strArr[i].length() > 0) {
                            try {
                                tmpV = Double.parseDouble(strArr[i]);
                            } catch (NumberFormatException e) {
                                continue;
                            }
//                System.out.println("\t\t"+tmpV);
                            if (tmpV != lastVal || (tmpV == lastVal && ++same < SAME_LIMIT)) {
                                if (tmpV != lastVal) same = 0;
                                lastVal = tmpV;
                                tmp++;
                                writer.write(tmpV + "\n");
                                if (tmp >= N) {
                                    writer.close();
                                    fos.close();
                                    txt_id++;
                                    tmp = 0;
                                    oFile = new File(folder + "\\N" + N + "\\tmp_" + txt_id + ".txt");
                                    fos = new BufferedOutputStream(new FileOutputStream(oFile));
                                    writer = new BufferedWriter(new OutputStreamWriter(fos, StandardCharsets.UTF_8), 5 * 1024 * 1024);
                                    writer.write("from: " + folder + "\n");
                                    break;
                                }
                            }
                        }
//            if(id>=50)break;
                }
                reader.close();
                fis.close();
                System.out.println("?\t\t" + tmp);
            }

    }


    public static void expand() throws IOException {
        int txt_id = 4;
        String filename = "liantong_data_from2018-12-19to2019-01-31_8260" + ".csv";
        String folder = "D:\\Study\\Lab\\iotdb\\real-world data\\云智慧";
        String filepath = folder + "\\" + filename;
        File file = new File(filepath);
        BufferedInputStream fis = null;
        fis = new BufferedInputStream(new FileInputStream(file));
        BufferedReader reader = new BufferedReader(new InputStreamReader(fis, StandardCharsets.UTF_8), 5 * 1024 * 1024);// 用5M的缓冲读取文本文件

        File oFile = new File(folder + "\\tmp_" + txt_id + ".txt");
        BufferedOutputStream fos = new BufferedOutputStream(new FileOutputStream(oFile));
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(fos, StandardCharsets.UTF_8), 5 * 1024 * 1024);
        writer.write("from: " + folder + " " + filename + "\n");

        String line = "";
        int id = 0;
        int tmp = 0;
        line = reader.readLine();//ignore first line
        double cycle[] = new double[1000000];
        while ((line = reader.readLine()) != null) {
            //TODO: write your business
            id++;
//            System.out.println("\t\t"+line);
            String[] strArr = line.split(",");
            double tmpV;
            try {
                tmpV = Double.parseDouble(strArr[1]);
            } catch (NumberFormatException e) {
                continue;
            }
            cycle[tmp++] = tmpV;
        }
        Random random = new Random();
        int tot = 0;
        for (int i = 0; i < N; i++)
            writer.write((cycle[i % tmp] + random.nextGaussian()) + "\n");

        reader.close();
        fis.close();
        System.out.println("?\t\t" + tmp);
        writer.close();
        fos.close();
    }

    public static void taxipredition8M() throws IOException {
        String dataset_name = "taxipredition8M";
        int txt_id = 4;
        String filename = "trainMehrak" + ".csv";
        String folder = "E:\\real-world data\\Kaggle\\taxipredictiondata8M";
        String filepath = folder + "\\" + filename;
        File file = new File(filepath);
        for (int i = 0; i < 4; i++) {
            File oFile = new File(folder + "\\" + dataset_name + "_" + i + ".txt");
            BufferedInputStream fis = null;
            fis = new BufferedInputStream(new FileInputStream(file));
            BufferedReader reader = new BufferedReader(new InputStreamReader(fis, StandardCharsets.UTF_8), 50 * 1024 * 1024);// 用50M的缓冲读取文本文件

            BufferedOutputStream fos = new BufferedOutputStream(new FileOutputStream(oFile));
            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(fos, StandardCharsets.UTF_8), 50 * 1024 * 1024);
            writer.write(dataset_name + " from: " + folder + " " + filename + "\n");

            String line = "";
            line = reader.readLine();//ignore first line
            while ((line = reader.readLine()) != null) {
//                System.out.println("\t\t"+line);
                String[] strArr = line.split(",");
                double tmpV;
                try {
                    tmpV = Double.parseDouble(strArr[3 + i]);
                    writer.write(tmpV + "\n");
                } catch (NumberFormatException e) {
                    continue;
                }
            }
            reader.close();
            fis.close();
            writer.close();
            fos.close();
        }
    }

    public static void physiological_stress_dataset() throws IOException {
        String dataset_name = "physiological_stress_dataset";
//        String folder = "E:\\real-world data\\Kaggle\\physiological stress dataset\\dataset\\0. interim\\swell\\eda";
        String folder = "E:\\real-world data\\Kaggle\\physiological stress dataset\\dataset\\0. interim\\wesad\\eda";
        File folderFile = new File(folder);
        File[] files = folderFile.listFiles();
        assert files != null;
        int txt_id = 0;


        File oFile = new File("E:\\real-world data\\Kaggle\\physiological stress dataset\\" + dataset_name + "_1.txt");
        BufferedOutputStream fos = new BufferedOutputStream(new FileOutputStream(oFile));
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(fos, StandardCharsets.UTF_8), 50 * 1024 * 1024);
        writer.write(dataset_name + " from: " + folder + "\n");

        for (File file : files)
            if (!file.isDirectory()) {
                System.out.println("\t\t\t" + file.getName());
                if (!file.getName().contains(".csv") && !file.getName().contains(".txt")) continue;
                BufferedInputStream fis = null;
                fis = new BufferedInputStream(new FileInputStream(file));
                BufferedReader reader = new BufferedReader(new InputStreamReader(fis, StandardCharsets.UTF_8), 50 * 1024 * 1024);

                String line = "";
                line = reader.readLine();//ignore first line
                while ((line = reader.readLine()) != null) {
                    String[] strArr = line.split(",");
                    double tmpV;
                    for (int i = 0; i < strArr.length; i++)
                        if (strArr[i].length() > 0) {
                            try {
                                tmpV = Double.parseDouble(strArr[i]);
                                writer.write(tmpV + "\n");
                            } catch (NumberFormatException e) {
                                continue;
                            }
                        }
                }
                reader.close();
                fis.close();
            }
        writer.close();
        fos.close();

    }


    public static void General_Practice_Prescribing_Data() throws IOException {
        String dataset_name = "General_Practice_Prescribing_Data";
        int txt_id = 4;
        String folder = "E:\\real-world data\\Kaggle\\General Practice Prescribing Data";
        File folderFile = new File(folder);
        File[] files = folderFile.listFiles();
        assert files != null;
        for (int i = 0; i < 2; i++) {
            File oFile = new File(folder + "\\" + dataset_name + "_" + i + ".txt");
            BufferedOutputStream fos = new BufferedOutputStream(new FileOutputStream(oFile));
            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(fos, StandardCharsets.UTF_8), 50 * 1024 * 1024);
            writer.write(dataset_name + " from: " + folder + "\n");
            for (File file : files)
                if (!file.isDirectory()) {
                    if (!file.getName().contains("+BNFT")) continue;
                    BufferedInputStream fis = null;
                    fis = new BufferedInputStream(new FileInputStream(file));
                    BufferedReader reader = new BufferedReader(new InputStreamReader(fis, StandardCharsets.UTF_8), 50 * 1024 * 1024);// 用50M的缓冲读取文本文件
                    String line = "";
                    line = reader.readLine();//ignore first line
                    while ((line = reader.readLine()) != null) {
//                System.out.println("\t\t"+line);
                        String[] strArr = line.split(",");
                        double tmpV;
                        try {
                            tmpV = Double.parseDouble(strArr[4 + i]);
                            writer.write(tmpV + "\n");
                        } catch (NumberFormatException e) {
                            continue;
                        }
                    }
                    reader.close();
                    fis.close();
//                    break;
                }
            writer.close();
            fos.close();
        }
    }


    public static void Spacecraft_Thruster_Firing_Tests_Dataset() throws IOException {
        String dataset_name = "SpacecraftThruster";
        int txt_id = 4;
        String folder = "D:\\Study\\Lab\\iotdb\\add_quantile_to_aggregation\\test_project_2\\Spacecraft Thruster Firing Tests Dataset\\dataset\\dataset\\train";
        File folderFile = new File(folder);
        File[] files = folderFile.listFiles();
        assert files != null;
        for (int i = 0; i < 1; i++) {
            File oFile = new File("./" + dataset_name + ".txt");
            BufferedOutputStream fos = new BufferedOutputStream(new FileOutputStream(oFile));
            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(fos, StandardCharsets.UTF_8), 50 * 1024 * 1024);
            writer.write(dataset_name + " from: " + "\"Spacecraft Thruster Firing Tests Dataset\" train_thrust" + "\n");
            for (File file : files)
                if (!file.isDirectory()) {
                    BufferedInputStream fis = null;
                    fis = new BufferedInputStream(new FileInputStream(file));
                    BufferedReader reader = new BufferedReader(new InputStreamReader(fis, StandardCharsets.UTF_8), 50 * 1024 * 1024);// 用50M的缓冲读取文本文件
                    String line = "";
                    line = reader.readLine();//ignore first line
                    while ((line = reader.readLine()) != null) {
//                System.out.println("\t\t"+line);
                        String[] strArr = line.split(",");
                        double tmpV;
                        try {
                            tmpV = Double.parseDouble(strArr[2 + i]);
                            writer.write(tmpV + "\n");
                        } catch (NumberFormatException e) {
                            continue;
                        }
                    }
                    reader.close();
                    fis.close();
//                    break;
                }
            writer.close();
            fos.close();
        }
    }
}