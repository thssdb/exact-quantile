import org.apache.datasketches.kll.KllDoublesSketch;
import org.apache.datasketches.kll.KllDoublesSketchIterator;

import java.io.*;
import java.util.Date;
import java.util.Random;

public class MainForMyKLL {
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

    public void testKLL_PAGE(){
        final int NNN=100000000;
        System.out.println("\t\t size:"+
            KllDoublesSketch.getMaxSerializedSizeBytes(45000, NNN, false));
        System.out.println("\t\t size:"+
            KllDoublesSketch.getMaxSerializedSizeBytes(45000, NNN/10, false));
        System.out.println("\t\t size:"+
            KllDoublesSketch.getMaxSerializedSizeBytes(10000, NNN, false));
        KllDoublesSketch worker = KllDoublesSketch.newHeapInstance(55);
        for(int i=0;i<NNN;i++)worker.update(i);
//        System.out.println(worker.toString(true, true));
//        worker.update(1e10);
//        System.out.println(worker.toString(true, false));
//        for(int i=0;i<14;i++)worker.update(i);
//        System.out.println(worker.toString(true, false));
//        for(int i=0;i<7;i++)worker.update(i);
//        System.out.println(worker.toString(true, true));
    }
    public void showKLLLevelSize(KllDoublesSketch worker){
        KllDoublesSketchIterator iterator = worker.iterator();
        int lastLevel=0,lastWeight=1,lastCnt=0;
        while(iterator.next()){
            if(iterator.getWeight()==lastWeight){
                lastCnt++;
            }else {
                while(iterator.getWeight()>lastWeight){
                    System.out.print("\t["+lastCnt+"]\t");
                    lastCnt=0;
                    lastWeight<<=1;
                }
                lastCnt=1;
            }
        }
        System.out.print("\t["+lastCnt+"]\n");
    }


    public void testMYKLL() throws IOException{
        final int maxMemoryByte=1<<20;
        final int maxSeriesByte=2048;
        final int NNN=275;
//        LongKLLSketch worker_1 = new LongKLLSketch(NNN, LEVEL_NUM, maxMemoryByte, maxSeriesByte);
//        for(int i=0;i<NNN;i++)worker_1.update(i);
////        worker_1.show();
//        worker_1.compactBeforeSerialization();
//        worker_1.show();

        long[] a = new long[NNN];
        for (int i = 0; i < NNN; i++) a[i] = i;
        for (int i = 1; i < NNN; i++){int p = random.nextInt(i);long tmp = a[i];a[i]=a[p];a[p]=tmp;}
        LongKLLSketch worker = new LongKLLSketch(NNN, maxMemoryByte, maxSeriesByte);
        for(int i=0;i<NNN;i++)worker.update(a[i]);
//        worker_2.show();
        worker.compactBeforeSerialization();
        worker.show();
        File f = new File("try_KLL.txt");
        OutputStream of = new FileOutputStream(f);
        worker.serialize(of);
        of.flush();

        InputStream i_f = new FileInputStream(f);

        worker = new LongKLLSketch(i_f);
        worker.show();
        worker.showLevelMaxSize();
        worker.showNum();

//        of = new FileOutputStream(f);
//        worker.serialize(of);
//        of.flush();
//        i_f = new FileInputStream(f);
//        worker = new LongKLLSketch(i_f);
//        worker.show();
//        worker.showLevelMaxSize();
//        worker.showNum();
//        System.out.println("\t\t"+worker.maxLevel);
//        System.out.println("\t\t ERR: "+calcMYDynamicKLLAvgMaxError(worker_2, NNN));
    }
    public double calcKLLAvgMaxError(KllDoublesSketch worker, int n){
        double avg = 0.0,mx=0;
        int[] bucket = new int[n];
        KllDoublesSketchIterator iterator = worker.iterator();
        while(iterator.next())
            bucket[(int)iterator.getValue()]+=iterator.getWeight();
        for(int i=1;i<n;i++)
            bucket[i]+=bucket[i-1];
        for(int i=1;i<n;i++){
            double delta = Math.abs(bucket[i-1]-(i));
            avg += delta/n/n;
            mx = Math.max(mx,delta);
        }
//        System.out.println("\t\t"+avg+"\t"+mx);
//        avg=mx=0;
//        for(int i=0;i<n;i++){
//            double delta = Math.abs(worker.getRank(a[i])*n-(i));
//            avg += delta/n/n;
//            mx = Math.max(mx,delta);
//        }System.out.println("\t\t"+avg+"\t"+mx);
        return mx+avg;
    }
//    public double calcMYKLLAvgMaxError(OldLongKLLSketch worker,int n){
//        double avg = 0.0,mx=0;
//        int[] bucket = new int[n];
//        int[] pos=worker.levelPos;long[] num = worker.num;
//        for(int level=0;level<worker.maxLevel;level++)
//            for(int i=pos[level];i<pos[level+1];i++)
//                bucket[(int)num[i]]+=1<<level;
//        for(int i=1;i<n;i++)
//            bucket[i]+=bucket[i-1];
//        for(int i=1;i<n;i++){
//            double delta = Math.abs(bucket[i-1]-(i));
//            avg += delta/n/n;
//            mx = Math.max(mx,delta);
//        }
//        return mx+avg;
//    }
    public double calcMYKLLAvgMaxError(LongKLLSketch worker,int n){
        double avg = 0.0,mx=0;
        int[] bucket = new int[n];
        int[] pos=worker.levelPos;long[] num = worker.num;
        for(int level=0;level<worker.cntLevel;level++)
            for(int i=pos[level];i<pos[level+1];i++)
                bucket[(int)num[i]]+=1<<level;
        for(int i=1;i<n;i++)
            bucket[i]+=bucket[i-1];
        for(int i=1;i<n;i++){
            double delta = Math.abs(bucket[i-1]-(i));
            avg += delta/n/n;
            mx = Math.max(mx,delta);
        }
    //        System.out.println("\t\t"+avg+"\t"+mx);
    //        avg=mx=0;
    //        for(int i=0;i<n;i++){
    //            double delta = Math.abs(worker.getApproxRank(a[i])-(i));
    //            avg += delta/n/n;
    //            mx = Math.max(mx,delta);
    //        }System.out.println("\t\t"+avg+"\t"+mx);
        return mx+avg;
    }
    public double calcMYHeapKLLAvgMaxError(HeapLongKLLSketch worker,int n){
        double avg = 0.0,mx=0;
        int[] bucket = new int[n];
        int[] pos=worker.levelPos;long[] num = worker.num;
        for(int level=0;level<worker.cntLevel;level++)
            for(int i=pos[level];i<pos[level+1];i++)
                bucket[(int)num[i]]+=1<<level;
        for(int i=1;i<n;i++)
            bucket[i]+=bucket[i-1];
        for(int i=1;i<n;i++){
            double delta = Math.abs(bucket[i-1]-(i));
            avg += delta/n/n;
            mx = Math.max(mx,delta);
        }
        return mx+avg;
    }

    public void testErrorCompare(){
        final int maxMemoryByte=1<<20;
        final int maxSeriesByte=1024*13;
        final int NNN=1600000;
        final int LEVEL_NUM=(int)Math.ceil(Math.log(NNN/(maxSeriesByte/8.0))/Math.log(2))+1;
//        LongKLLSketch worker3 = new LongKLLSketch(NNN_CHUNK, 11);
        long[] a = new long[NNN];
        for (int i = 0; i < NNN; i++) a[i] = i;
        double err1=0,err1_SERI=0,err2=0,err2_SERI=0,avg1=0,avg1_SERI=0,avg2=0,avg2_SERI=0,mx1=0,mx1_SERI=0,mx2=0,mx2_SERI=0;
        int TEST_CASE = (int)(6e7/NNN);Random random = new Random();
        int KLL_SERI_K = 0;
        for(int t=0;t<TEST_CASE;t++) {
            for (int i = 1; i < NNN; i++)/*if((i&15)<=12)*/{
                int p = random.nextInt(i);
                long tmp = a[i];a[i]=a[p];a[p]=tmp;
            }
//            if(t==0)for(int i=0;i<NNN;i++)a[i]=i;
//            for(int i=0;i<10;i++)System.out.print("\t"+a[i]);System.out.println();

            KllDoublesSketch worker1;int KLL_K=0;
            for(int k=1<<20;k>0;k>>>=1)if(KllDoublesSketch.getMaxSerializedSizeBytes(KLL_K+k, NNN, false)<=maxMemoryByte)
                KLL_K+=k;
            KLL_K = Math.min(KLL_K,NNN);
            worker1=KllDoublesSketch.newHeapInstance(KLL_K);
            LongKLLSketch worker2 = new LongKLLSketch(NNN, maxMemoryByte, maxSeriesByte);
            for (int i = 0; i < NNN; i++)worker1.update(a[i]);
            for (int i = 0; i < NNN; i++){
                worker2.update(a[i]);
            }
//            Arrays.sort(a, 0, NNN);
//            System.out.println(worker1.toString(true, true));

//            worker2.show();
//        System.out.println("\t\t\t?? "+worker2.getApproxRank(3500));

            err1 = calcKLLAvgMaxError(worker1, NNN);
            avg1 = avg1+(err1-Math.floor(err1))/TEST_CASE;
            mx1 = Math.max(mx1, Math.floor(err1));

            KllDoublesSketch worker1_SERI;
            if(KLL_SERI_K==0) {
                for(int i=0;i<2;i++)
                for (int k = 1 << Math.min(15, (int) Math.floor(Math.log(maxSeriesByte / 8.0) / Math.log(2))); k > 0; k >>>= 1) {
                    if (KLL_SERI_K == 0 && k < 8) break;
                    worker1_SERI = KllDoublesSketch.newHeapInstance(KLL_SERI_K + k);
                    worker1_SERI.merge(worker1);
                    if (worker1_SERI.getSerializedSizeBytes() <= maxSeriesByte) KLL_SERI_K += k;
                }
                if (KLL_SERI_K == 0) KLL_SERI_K = 8;
            }
            worker1_SERI = KllDoublesSketch.newHeapInstance(KLL_SERI_K);
            worker1_SERI.merge(worker1);

            err1_SERI = calcKLLAvgMaxError(worker1_SERI, NNN);
            avg1_SERI = avg1_SERI+(err1_SERI-Math.floor(err1_SERI))/TEST_CASE;
            mx1_SERI = Math.max(mx1_SERI, Math.floor(err1_SERI));

            err2 = calcMYKLLAvgMaxError(worker2, NNN);
            avg2 = avg2+(err2-Math.floor(err2))/TEST_CASE;
            mx2 = Math.max(mx2, Math.floor(err2));
            if(t+1==TEST_CASE){
//                System.out.println(worker1.toString(true, false));
//                System.out.println(worker1_SERI.toString(true, false));
                showKLLLevelSize(worker1);
                showKLLLevelSize(worker1_SERI);
                worker2.show();
            }

            worker2.compactBeforeSerialization();
            err2_SERI = calcMYKLLAvgMaxError(worker2, NNN);
            avg2_SERI = avg2_SERI+(err2_SERI-Math.floor(err2_SERI))/TEST_CASE;
            mx2_SERI = Math.max(mx2_SERI, Math.floor(err2_SERI));
//            System.out.println(worker1.toString(true, true));
//            worker2.show();
//            for(int i=0;i<NNN;i+=NNN/10)System.out.println("\t\t[ERROR] "+worker2.getApproxRank(i));
            if(t+1==TEST_CASE){
                worker2.show();
            }
            if((t*10)%TEST_CASE<10)
                System.out.println("\t\t[MIDDLE] KLL:" + "avg:"+avg1+" mx:"+mx1 + "\t\t" + "KLL_SERI:" + "avg:"+avg1_SERI+" mx:"+mx1_SERI + "\t||\t" + "myKLL:" + "avg:"+avg2+" mx:"+mx2 + "\t\t" + "myKLL_SERI:" + "avg:"+avg2_SERI+" mx:"+mx2_SERI);
        }avg1*=NNN;avg1_SERI*=NNN;avg2*=NNN;avg2_SERI*=NNN;
        System.out.println("\t\t\t[parameter]\tmaxMem:"+maxMemoryByte+"\tmaxSERI:"+maxSeriesByte);
        System.out.println("\t\t[ERROR] KLL:" + "avg:"+avg1+" mx:"+mx1 + "\t\t" + "KLL_SERI:" + "avg:"+avg1_SERI+" mx:"+mx1_SERI + "\t||\t" + "myKLL:" + "avg:"+avg2+" mx:"+mx2 + "\t\t" + "myKLL_SERI:" + "avg:"+avg2_SERI+" mx:"+mx2_SERI);

    }
    public void testSeriDeseri() throws IOException {
        final int maxMemoryByte=1<<20;
        final int maxSeriesByte=1024*13;
        final int NNN=1600000;

        LongKLLSketch worker = new LongKLLSketch(NNN, maxMemoryByte, maxSeriesByte);
        for(int i=0;i<NNN/3;i++)worker.update(i);
        File f = new File("try_seri.txt");
        OutputStream of = new FileOutputStream(f);
        worker.serialize(of);
        of.flush();
        of.close();
        InputStream inf = new FileInputStream(f);
        LongKLLSketch worker2 = new LongKLLSketch(inf, maxMemoryByte, maxSeriesByte);
        inf.close();
//        worker.show();
//        worker.showLevelMaxSize();
//        worker2.show();
//        worker2.showLevelMaxSize();
        for(int i=0;i<NNN/3;i++)worker2.update(i);
        for(int i=0;i<NNN/3;i++)worker.update(i);
        worker.show();
//        worker.showLevelMaxSize();
        worker2.show();
//        worker2.showLevelMaxSize();
    }
    public void testMerge() throws IOException {
        final int maxMemoryByte=1<<20;
        final int maxSeriesByte=480/*1024*13*/;
        final int NNN=7000;

        LongKLLSketch worker = new LongKLLSketch(NNN, maxMemoryByte, maxSeriesByte);
        for(int i=0;i<NNN;i++)worker.update(i);
        File f = new File("try_seri.txt");
        OutputStream of = new FileOutputStream(f);
        worker.serialize(of);
        of.flush();
        of.close();
        InputStream inf = new FileInputStream(f);
        LongKLLSketch worker2 = new LongKLLSketch(inf, maxMemoryByte, maxSeriesByte);
        inf.close();

        HeapLongKLLSketch worker_3 = new HeapLongKLLSketch(1<<20);
        for(int i=0;i<10000;i++)worker_3.merge(worker2);
        worker_3.show();
        worker_3.showLevelMaxSize();
        for(int i=0;i<10;i++)System.out.println("\t\t"+worker_3.getApproxRank(350+i*700));
    }

    public void testMergeErr(final int TEST_CASE, final int MERGE_NUM){
        final int maxMemoryByte=/*1<<20*/480*100+233;
        final int maxSeriesByte=565;
        final int NNN=7989;
//        LongKLLSketch worker3 = new LongKLLSketch(NNN_CHUNK, 11);
        long[] a = new long[NNN];
//        double err1=0,err1_SERI=0,err2=0,err2_SERI=0,avg1=0,avg1_SERI=0,avg2=0,avg2_SERI=0,mx1=0,mx1_SERI=0,mx2=0,mx2_SERI=0;
        double merge_err1=0,merge_err2=0,merge_avg1=0,merge_avg2=0,merge_mx1=0,merge_mx2=0;
        final int MERGE_NNN = NNN*MERGE_NUM;Random random = new Random();
        long[] merge_a = new long[MERGE_NNN];
        for (int i = 0; i < MERGE_NNN; i++) merge_a[i] = i;

        for(int t=0;t<TEST_CASE;t++) {
            for (int i = 1; i < MERGE_NNN; i++){
                int p = random.nextInt(i);
                long tmp = merge_a[i];merge_a[i]=merge_a[p];merge_a[p]=tmp;
            }
            int KLL_SERI_K = 0;
            int MERGE_KLL_K=0;KllDoublesSketch mergeWorker1=null;
            HeapLongKLLSketch mergeWorker2 = new HeapLongKLLSketch(maxMemoryByte);
            for (int num = 0; num < MERGE_NUM; num++) {

//                KllDoublesSketch worker1;
                int KLL_K = 0;
                for (int k = 1 << 20; k > 0; k >>>= 1)
                    if(KLL_K + k<65536)
                    if (KllDoublesSketch.getMaxSerializedSizeBytes(KLL_K + k, NNN, false) <= maxMemoryByte)
                        KLL_K += k;
                KLL_K = Math.min(KLL_K, NNN);
//                worker1 = KllDoublesSketch.newHeapInstance(KLL_K);
                LongKLLSketch worker2 = new LongKLLSketch(NNN, maxMemoryByte, maxSeriesByte);
                for (int i = 0; i < NNN; i++) a[i] = merge_a[NNN * num + i];
//                for (int i = 0; i < NNN; i++) worker1.update(a[i]);
                for (int i = 0; i < NNN; i++) {
                    worker2.update(a[i]);
                }

//                KllDoublesSketch worker1_SERI;
//                if (KLL_SERI_K == 0) {
//                    for (int k = 1 << Math.min(15, (int) Math.floor(Math.log(maxSeriesByte / 8.0) / Math.log(2))); k > 0; k >>>= 1) {
//                        if (KLL_SERI_K == 0 && k < 8) break;
//                        worker1_SERI = KllDoublesSketch.newHeapInstance(KLL_SERI_K + k);
//                        worker1_SERI.merge(worker1);
//                        if (worker1_SERI.getSerializedSizeBytes() <= maxSeriesByte) KLL_SERI_K += k;
//                    }
//                    worker1_SERI = KllDoublesSketch.newHeapInstance(KLL_SERI_K);
//                    worker1_SERI.merge(worker1);
//                    for (int k = 1 << 20; k > 0; k >>>= 1)
//                        if (MERGE_KLL_K + k < 65536) {
//                            mergeWorker1 = KllDoublesSketch.newHeapInstance(MERGE_KLL_K + k);
//                            boolean flag=true;
//                            for (int tt = 0; tt < MERGE_NUM; tt++) {
//                                mergeWorker1.merge(worker1_SERI);
//                                flag&=mergeWorker1.getSerializedSizeBytes() <= maxMemoryByte;
//                            }
////                            System.out.println("\t\t??\t\t"+(MERGE_KLL_K + k)+":"+mergeWorker1.getSerializedSizeBytes());
//                            if (flag) MERGE_KLL_K += k;
//                        }
//                }
//                if (num == 0)
//                    mergeWorker1 = KllDoublesSketch.newHeapInstance(MERGE_KLL_K);
//                worker1_SERI = KllDoublesSketch.newHeapInstance(KLL_SERI_K);
//                worker1_SERI.merge(worker1);

//                mergeWorker1.merge(worker1_SERI);

                worker2.compactBeforeSerialization();
                mergeWorker2.mergeWithTempSpace(worker2);
                if(t+1==TEST_CASE&&num+1==MERGE_NUM) {
//                    showKLLLevelSize(worker1_SERI);
                    worker2.show();
                    worker2.showLevelMaxSize();
                    worker2.showNum();
//                    System.out.println("\t\t\t[parameter]\tmaxMem:" + maxMemoryByte + "\tmaxSERI:" + maxSeriesByte);
//                    showKLLLevelSize(mergeWorker1);
//                    mergeWorker2.show();
//                    System.out.println("\t\t\t\t"+mergeWorker1.toString(true,false));
                }
            }
//            merge_err1 = calcKLLAvgMaxError(mergeWorker1, MERGE_NNN);
//            merge_avg1 = merge_avg1 + (merge_err1 - Math.floor(merge_err1)) / TEST_CASE;
//            merge_mx1 = Math.max(merge_mx1, Math.floor(merge_err1));

            merge_err2 = calcMYHeapKLLAvgMaxError(mergeWorker2, MERGE_NNN);
            merge_avg2 = merge_avg2 + (merge_err2 - Math.floor(merge_err2)) / TEST_CASE;
            merge_mx2 = Math.max(merge_mx2, Math.floor(merge_err2));
        }


//        System.out.println("\t\t[ERROR] MERGE_KLL:" + "avg:"+merge_avg1*MERGE_NNN+" mx:"+merge_mx1 + "\t||\t" + "MERGE_MYKLL:" + "avg:"+merge_avg2*MERGE_NNN+" mx:"+merge_mx2);
//        System.out.print("\t("+Math.ceil(merge_avg2*MERGE_NNN)+","+Math.ceil(merge_mx2)+")");
        System.out.print(MERGE_NUM+"\t"+Math.ceil(merge_avg2*MERGE_NNN)+"\t"+Math.ceil(merge_mx2)+"\n");
//        System.out.println("");
    }

    public void testMYHeapKLL(final int TEST_CASE, final int NNN){
        final int maxMemoryByte=480*100+233;
//        LongKLLSketch worker_1 = new LongKLLSketch(NNN, LEVEL_NUM, maxMemoryByte, maxSeriesByte);
//        for(int i=0;i<NNN;i++)worker_1.update(i);
////        worker_1.show();
//        worker_1.compactBeforeSerialization();
//        worker_1.show();
        double merge_err2,merge_avg2=0,merge_mx2=0;

        long[] a = new long[NNN];
        for (int i = 0; i < NNN; i++) a[i] = i;
        for(int t=0;t<TEST_CASE;t++) {
            for (int i = 1; i < NNN; i++) {
                int p = random.nextInt(i);
                long tmp = a[i];
                a[i] = a[p];
                a[p] = tmp;
            }
            HeapLongKLLSketch worker = new HeapLongKLLSketch(maxMemoryByte);
            for (int i = 0; i < NNN; i++) worker.update(a[i]);
            merge_err2 = calcMYHeapKLLAvgMaxError(worker, NNN);
            merge_avg2 = merge_avg2 + (merge_err2 - Math.floor(merge_err2)) / TEST_CASE;
            merge_mx2 = Math.max(merge_mx2, Math.floor(merge_err2));
            if(t+1==TEST_CASE){
//                worker.show();
            }
        }
//        System.out.print("\t("+merge_avg2*NNN+","+Math.ceil(merge_mx2)+")");
        System.out.print(NNN/7000+"\t"+merge_avg2*NNN+"\n");
    }


    public static void testMyKLLRank(final int NNN){
        LongKLLSketch worker = new LongKLLSketch(NNN,1<<20,560);
        for(int i=0;i<NNN;i++)worker.update(i);
        worker.compactBeforeSerialization();
        worker.show();
        worker.showLevelMaxSize();
        worker.showNum();
    }



    public void testMYStrictKLL() throws IOException{
        final int maxMemoryByte=1<<12;
        final int NNN=1<<24;
        long[] a = new long[NNN];
        for (int i = 0; i < NNN; i++) a[i] = i;
        for (int i = 1; i < NNN; i++){int p = random.nextInt(i);long tmp = a[i];a[i]=a[p];a[p]=tmp;}


        long TIME_2 = new Date().getTime();
        KLLSketchForQuantile worker_strict = new HeapLongStrictKLLSketch(maxMemoryByte);
        for(int i=0;i<NNN;i++)worker_strict.update(a[i]);

        long TIME_1 = new Date().getTime();
        KLLSketchForQuantile worker = new HeapLongKLLSketch(maxMemoryByte);
        for(int i=0;i<NNN;i++)worker.update(a[i]);
        TIME_1 = new Date().getTime()-TIME_1;


        TIME_2 = new Date().getTime()-TIME_2;
//        worker.show();
        worker.show();
        worker_strict.show();
        double avgE=0,avgE_2 =0;
        for(double q=0.01;q<=0.991;q+=0.01){
            long rk = worker.getApproxRank((int)(q*NNN));
            long rk_2 = worker_strict.getApproxRank((int)(q*NNN));
            System.out.print("\t\t"+rk+
                "\t\t"+rk_2);
            double err = 1.0*(rk-(int)(q*NNN))/NNN;
            double err_2 = 1.0*(rk_2-(int)(q*NNN))/NNN;
            System.out.println("\t\t"+err+
                "\t\t"+err_2);
            avgE+=err/99;
            avgE_2+=err_2/99;
        }
        System.out.println("\t\t\t"+TIME_1+"\t"+TIME_2);
        System.out.println("\t\t\t"+avgE+"\t"+avgE_2);
    }


    public static void main(String[] args) throws IOException{
        MainForMyKLL main;
        main = new MainForMyKLL();
        main.testMYStrictKLL();
//        main.testKLL_PAGE();System.out.println("\n\n");
//       for(int i=0;i<100;i++) main.testMYKLL();System.out.println("\n\n");
//        main.testErrorCompare();
//        main.testSeriDeseri();
//        main.testMerge();
//        main.testMergeErr(1,1);
//        for(int num=1;num<=9;num+=1){
//            main.testMergeErr(16,num);
//        }
//        for(int num=10;num<=100;num+=10){
//            main.testMergeErr(100,num);
//        }
//        for(int num=2700;num<=2700;num+=500){
//            main.testMergeErr(200,num);
//        }
//        for(int test=1<<8;test<=1<<13;test*=2){
//            main.testMergeErr(test,100);
//        }
//        for(int i=1;i<=10;i++)
//            main.testMYHeapKLL(3000, i*7000);
//        for(int i=20;i<=100;i+=10)
//            main.testMYHeapKLL(500, i*7000);
//        for(int i=200;i<=2000;i+=200)
//            main.testMYHeapKLL(60, i*7000);
//        testMyKLLRank(8704);
    }
}

// 	[8478]		[1]		[16986]		[55203]		[0]		[33811]
//	(14.29260147472178,81.0)
//	[7872]		[0]		[17326]		[13445]		[88454]
//	(12.564914867994048,76.0)