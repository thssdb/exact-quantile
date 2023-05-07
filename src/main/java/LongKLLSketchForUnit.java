import com.google.common.primitives.UnsignedBytes;
import org.apache.iotdb.tsfile.utils.ReadWriteIOUtils;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;

// based on KLL Sketch in DataSketch. See https://github.com/apache/datasketches-java/tree/master/src/main/java/org/apache/datasketches/kll
// This is an implementation for long data type.
// We changed the behaviour of serialization to reduce the serialization size. In our situation, there will be few update operations after deserialization.
public class LongKLLSketchForUnit extends KLLSketchForQuantile{
  public LongKLLSketchForUnit(int N, int cntLevel,int sketchSize) {
    this.cntLevel = cntLevel;
    this.N = N;
    levelPos = new int[cntLevel+1];
    num = new long[sketchSize];
  }
  public void addWhenDivide(int level,long v){
    num[levelPos[level+1]++]=v;
  }


  public String toString(){
    final StringBuilder sb = new StringBuilder();
    sb.append(N);
    sb.append((byte)cntLevel);
    sb.append((short)(levelPos[cntLevel]-levelPos[0]));
    for(int i=0;i<cntLevel;i++)sb.append(levelPos[i]);
    return sb.toString();
  }

  public int serialize(OutputStream outputStream) throws IOException {
    int byteLen = 0;
    byteLen+= ReadWriteIOUtils.write((int)N, outputStream);
    byteLen+=ReadWriteIOUtils.write((byte)cntLevel, outputStream);
    int numLEN = getNumLen();
    if(numLEN<256)
      for(int i=0;i<=cntLevel;i++)byteLen+=ReadWriteIOUtils.write((byte)(levelPos[i]), outputStream);
    else
      for(int i=0;i<=cntLevel;i++)byteLen+=ReadWriteIOUtils.write((short)(levelPos[i]), outputStream);
    for(int i=levelPos[0];i<levelPos[cntLevel];i++)byteLen+=ReadWriteIOUtils.write(num[i],outputStream);
    return byteLen;
  }
  public LongKLLSketchForUnit(InputStream inputStream) throws IOException {
    this.N = ReadWriteIOUtils.readInt(inputStream);
    this.cntLevel = ReadWriteIOUtils.readByte(inputStream);
    int numLEN = ReadWriteIOUtils.readShort(inputStream);
    num = new long[numLEN];
    level0Sorted = true;
    levelPos = new int[cntLevel];
    if (numLEN < 256)
      for (int i = 0; i <= cntLevel; i++)
        levelPos[i] = UnsignedBytes.toInt(ReadWriteIOUtils.readByte(inputStream));
    else
      for (int i = 0; i <= cntLevel; i++)
        levelPos[i] = ReadWriteIOUtils.readShort(inputStream);

    for (int i = 0; i < numLEN; i++) num[levelPos[0]+i] = ReadWriteIOUtils.readLong(inputStream);
  }

  public LongKLLSketchForUnit(ByteBuffer byteBuffer) {
    this.N = ReadWriteIOUtils.readInt(byteBuffer);
    this.cntLevel = ReadWriteIOUtils.readByte(byteBuffer);
    int numLEN = ReadWriteIOUtils.readShort(byteBuffer);
    num = new long[numLEN];
    level0Sorted = true;
    levelPos = new int[cntLevel];
    if (numLEN < 256)
      for (int i = 0; i <= cntLevel; i++)
        levelPos[i] = UnsignedBytes.toInt(ReadWriteIOUtils.readByte(byteBuffer));
    else
      for (int i = 0; i <= cntLevel; i++)
        levelPos[i] = ReadWriteIOUtils.readShort(byteBuffer);

    for (int i = 0; i < numLEN; i++) num[levelPos[0]+i] = ReadWriteIOUtils.readLong(byteBuffer);
  }
}