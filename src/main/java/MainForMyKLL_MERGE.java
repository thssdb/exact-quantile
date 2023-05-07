import java.io.*;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Random;

public class MainForMyKLL_MERGE {
    static int N = 20000000;
    final long KK = 47, num_mask = 65536;

    long seed;
    Random random = new Random(233);
    private long nextLong(){
        return random.nextLong();
//        seed ^= (seed << 21);
//        seed ^= (seed >>> 35);
//        seed ^= (seed << 4);
//        return seed;
    }
    private void reset(){seed=2333L;random.setSeed(233);}
    private long nextNumber(long x, int i) {
        long num=x;
        num=x&262143;
        if((i&15)<=3)
            num = nextLong();
        return num;
    }



    public void testMerge() throws IOException {
        final int maxMemoryByte=1<<19;
        final int maxSeriesByte=545/*1024*13*/;
        final int NNN=7989;

//        LongKLLSketch worker = new LongKLLSketch(NNN, maxMemoryByte, maxSeriesByte);
//        for(int i=0;i<NNN;i++)worker.update(i);
        File f = new File("try_seri.txt");
//        OutputStream of = new FileOutputStream(f);
//        worker.serialize(of);
//        of.flush();
//        of.close();
        InputStream inf = new FileInputStream(f);
        LongKLLSketch worker2 = new LongKLLSketch(inf, maxMemoryByte, maxSeriesByte);
        inf.close();

        worker2.show();
        final int merge_num = 12000,list_len = 1000,list_num = merge_num/list_len;
        long ST_TIME = new Date().getTime(),TIME_1,TIME_2;
        HeapLongKLLSketch worker_3 = new HeapLongKLLSketch(maxMemoryByte);
        List<KLLSketchForQuantile> workerList = new ArrayList<>();
        for(int i=0;i<list_len;i++)workerList.add(worker2);
        for(int i=0;i<list_num;i++)worker_3.mergeWithTempSpace(workerList);
        worker_3.show();
        worker_3.showLevelMaxSize();
        for(int i=0;i<10;i++)System.out.println("\t\t"+worker_3.getApproxRank(350+i*700));
        TIME_1 = new Date().getTime()-ST_TIME;

        ST_TIME = new Date().getTime();
        HeapLongKLLSketch worker_4 = new HeapLongKLLSketch(maxMemoryByte);
        for(int i=0;i<merge_num;i++)worker_4.mergeWithTempSpace(worker2);
        worker_4.show();
        worker_4.showLevelMaxSize();
        for(int i=0;i<10;i++)System.out.println("\t\t"+worker_4.getApproxRank(350+i*700));
        TIME_2 = new Date().getTime()-ST_TIME;
        System.out.println("TIME:\t\t\t"+TIME_1+"\t\t"+TIME_2);
    }




    public static void main(String[] args) throws IOException{
        MainForMyKLL_MERGE main;
        main = new MainForMyKLL_MERGE();
        main.testMerge();
    }
}

// 	[8478]		[1]		[16986]		[55203]		[0]		[33811]
//	(14.29260147472178,81.0)
//	[7872]		[0]		[17326]		[13445]		[88454]
//	(12.564914867994048,76.0)