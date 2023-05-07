import java.io.*;
import java.util.Date;
import java.util.Random;

public class MainForTDigestForPageStat {
    static int N = 20000000;
    final long KK = 47, num_mask = 65536;

    long seed;
    Random random = new Random(233);

    public void test_stat_writeAndRead(int stat_N,short stat_K) throws IOException{
        TDigestForPageStatWriter worker = new TDigestForPageStatWriter(stat_N,stat_K);
        for(int i=0;i<stat_N;i++)worker.update(/*(long)Math.pow(1.5,i)*/random.nextLong());

        File f = new File("try_seri_tdigest.txt");
        OutputStream of = new FileOutputStream(f);
        worker.serialize(of);
        of.flush();
        of.close();
        InputStream inF = new FileInputStream(f);
        TDigestForPageStatReader reader = new TDigestForPageStatReader(inF);
        inF.close();
    }

    public void test_TDigest_calc(int stat_N) throws IOException {
        TDigestForPageStatCalculator worker = new TDigestForPageStatCalculator(32768*4*8,0.5);
        long ST=new Date().getTime();
        long[] a = new long[stat_N];
        for(int i=0;i<stat_N;i++)a[i] = i;
//        for(int i=0;i<stat_N;i++)a[i] = random.nextLong();
        for (int i = 1; i < stat_N; i++){int p = random.nextInt(i);long tmp = a[i];a[i]=a[p];a[p]=tmp;}
        for(int i=0;i<stat_N;i++)
            worker.updateFromData(a[i]);
        System.out.println("\t\ttime:"+(new Date().getTime()-ST));

        worker.show();
        worker.showCluster();
        System.out.println("\t\t\t??\t\t"+worker.minValueWithRank(50)+" "+worker.maxValueWithRank(50));

//        Arrays.sort(a);
//        double err=0,err_c=0;
//        for(int i=0;i<stat_N;i+=stat_N/300){
//            long v1=worker.minValueWithRank(i),v2=worker.maxValueWithRank(i);
//            long ans=v1+((v2-v1)>>>1),rk=0;
//            for(int j=0;j<stat_N;j++)if(a[j]<=ans)rk++;
//            err+=Math.abs(rk-i);
//            err_c++;
//        }
//        System.out.println("\t\tavg_err:"+err/err_c);

//        File f = new File("try_seri_tdigest.txt");
//        InputStream inF = new FileInputStream(f);
//        TDigestForPageStatReader reader = new TDigestForPageStatReader(inF);
//        inF.close();
//        worker.updateFromPageStat(reader);
//        worker.show();
//        worker.showCluster();
    }

    public void test_stat_IOCalc(int stat_N,short stat_K) throws IOException{
        TDigestForPageStatWriter worker = new TDigestForPageStatWriter(stat_N,stat_K);
        long[] a = new long[stat_N];
        for(int i=0;i<stat_N;i++)a[i] = i;
        for (int i = 1; i < stat_N; i++){int p = random.nextInt(i);long tmp = a[i];a[i]=a[p];a[p]=tmp;}
        for(int i=0;i<stat_N;i++)worker.update(a[i]);

        File f = new File("try_seri_tdigest.txt");
        OutputStream of = new FileOutputStream(f);
        worker.serialize(of);
        of.flush();
        of.close();
        InputStream inF = new FileInputStream(f);
        TDigestForPageStatReader reader = new TDigestForPageStatReader(inF);
        inF.close();
        TDigestForPageStatCalculator calc = new TDigestForPageStatCalculator(32768*4*8,0.5);
        calc.updateFromPageStat(reader);
        calc.show();
    }

    public static void main(String[] args) throws IOException{
        MainForTDigestForPageStat main;
        main = new MainForTDigestForPageStat();
//        main.test_stat_writeAndRead(7989,(short)256);
        main.test_TDigest_calc(100);
//        main.test_stat_IOCalc(7989,(short)256);
    }
}