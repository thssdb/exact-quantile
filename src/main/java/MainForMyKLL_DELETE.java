import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class MainForMyKLL_DELETE {
    static int N = 20000000;
    final long KK = 47, num_mask = 65536;

    long seed;
    Random random = new Random(233);

    public void testDelete(int N_TOT,int N_STAT,double del_rate,int TEST_CASE) throws IOException {
        final int maxMemoryByte=1<<20;
        final int maxSeriesByte=549/*1024*13*/;


        boolean[] deleted = new boolean[N_TOT];
        long[] a = new long[N_TOT];
        for(int i=0;i<N_TOT;i++)a[i]=i;
        for(int i=1;i<N_TOT;i++){int p = random.nextInt(i);long tmp = a[i];a[i]=a[p];a[p]=tmp;}
        int[] value_rank = new int[N_TOT];
        for(int i=0;i<N_TOT;i++)deleted[i] = /*random.nextFloat()<del_rate*//*a[i]*/i<del_rate*N_TOT;
        for(int i=0;i<N_TOT;i++)value_rank[(int)a[i]]=!deleted[i]?1:0;
        for(int i=1;i<N_TOT;i++)value_rank[i]+=value_rank[i-1];

        File f = new File("try_delete.txt");

//        OutputStream of = new FileOutputStream(f);
//        LongKLLSketch worker = new LongKLLSketch(N_STAT, maxMemoryByte, maxSeriesByte);
//        for(int i=0;i<N_TOT;i++) {
//            worker.update(a[i]);
//            if(worker.getN()==N_STAT || i+1==N_TOT) {
//                worker.serialize(of);
//                worker = new LongKLLSketch(N_STAT, maxMemoryByte, maxSeriesByte);
//            }
//        }
//        of.flush();
//        of.close();
        InputStream inf = new FileInputStream(f);
        List<LongKLLSketch> stat_list = new ArrayList<>();
        int stat_num = (N_TOT-1)/N_STAT+1;
        for(int i=0;i<stat_num;i++) {
            LongKLLSketch tmpWorker = new LongKLLSketch(inf, maxMemoryByte, maxSeriesByte);
            stat_list.add(tmpWorker);
        }
        inf.close();

//        for(int i=0;i<5;i++)System.out.print(a[i]+"\t");System.out.println();

        double test_err_normal = 0, test_err_merge_delete = 0;
        for(int T=0;T<TEST_CASE;T++) {
            HeapLongKLLSketch KLL_NORMAL = new HeapLongKLLSketch(maxMemoryByte);
            for (int i = 0; i < N_TOT; i++) if (!deleted[i]) KLL_NORMAL.update(a[i]);

            int mem_for_merge = maxMemoryByte / 3;
            HeapLongKLLSketch KLL_MERGE = new HeapLongKLLSketch(mem_for_merge);
            int stat_mem = 560;
            int merge_num = mem_for_merge / stat_mem;
            for (int i = 0; i < stat_num; i += merge_num) {
                List<KLLSketchForQuantile> tmpList = new ArrayList<>();
                for (int j = i; j < i + merge_num && j < stat_num; j++)
                    tmpList.add(stat_list.get(j));
                KLL_MERGE.mergeWithTempSpace(tmpList);
            }
            HeapLongKLLSketch KLL_DELETED = new HeapLongKLLSketch(mem_for_merge);
            for (int i = 0; i < N_TOT; i++) if (deleted[i]) KLL_DELETED.update(a[i]);

//            KLL_NORMAL.show();
//            KLL_NORMAL.showLevelMaxSize();
//            KLL_MERGE.show();
//            KLL_MERGE.showLevelMaxSize();
//            KLL_DELETED.show();
//            KLL_DELETED.showLevelMaxSize();

            int query_interval = N_TOT / 1000, query_times = 0;
            long actual_N = KLL_NORMAL.getN();
            double err_normal = 0, err_merge_delete = 0;
            for (int query_x = random.nextInt(query_interval); query_x < N_TOT; query_x += query_interval) {
                query_times++;
                int actual_rank = value_rank[query_x] - 1;
                int normal_rank = KLL_NORMAL.getApproxRank(query_x);
                int merge_delete_rank = KLL_MERGE.getApproxRank(query_x) - KLL_DELETED.getApproxRank(query_x);
                err_normal += 1.0 * Math.abs(actual_rank - normal_rank) / actual_N;
                err_merge_delete += 1.0 * Math.abs(actual_rank - merge_delete_rank) / actual_N;
//                System.out.println("\t\t" + actual_rank + "\t" + normal_rank + "\t" + merge_delete_rank);
            }
            err_normal /= query_times;
            err_merge_delete /= query_times;
//            System.out.println("\t\t\t" + del_rate + "\t" + err_normal + "\t" + err_merge_delete);
            test_err_normal+=err_normal;test_err_merge_delete+=err_merge_delete;
        }
        System.out.println("\t\t\t" + del_rate + "\t" + test_err_normal/TEST_CASE + "\t" + test_err_merge_delete/TEST_CASE);
    }




    public static void main(String[] args) throws IOException{
        MainForMyKLL_DELETE main;
        main = new MainForMyKLL_DELETE();
        for(double rate=0;rate<=0.3;rate+=0.02)
            main.testDelete(100000000,7989, rate,16);
    }
}

// 	[8478]		[1]		[16986]		[55203]		[0]		[33811]
//	(14.29260147472178,81.0)
//	[7872]		[0]		[17326]		[13445]		[88454]
//	(12.564914867994048,76.0)