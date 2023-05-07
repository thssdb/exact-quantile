/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

import org.eclipse.collections.api.map.primitive.MutableLongLongMap;
import org.eclipse.collections.api.tuple.primitive.LongLongPair;
import org.eclipse.collections.impl.map.mutable.primitive.LongLongHashMap;

public class GSHashMapForStat {
  int bitsOfValue, remainingBits, maxSize;

  MutableLongLongMap hashMap;

  public GSHashMapForStat(int bits, int maxSize) {
    bitsOfValue = remainingBits = bits;
    hashMap = new LongLongHashMap(maxSize);
    this.maxSize = maxSize;
  }
  public GSHashMapForStat(int bits, int maxSize, int remaining) {
    this(bits,maxSize);
    remainingBits = remaining;
  }

  private void rebuild() {
//    System.out.println("[rebuild]+remaining:"+remainingBits+"  tot:"+totSize);
    int SHR = 1;
    MutableLongLongMap newMap = new LongLongHashMap(maxSize);
    hashMap.forEachKeyValue((k,v)->newMap.addToValue(k>>>1,v));

    if (newMap.size() >= maxSize) {
      SHR = 0;
      for(int delta = Integer.highestOneBit(remainingBits-3);delta>0;delta>>>=1)
      if(SHR+delta<=remainingBits-3){
        newMap.clear();
        final int shr = SHR+delta;
        hashMap.forEachKeyValue((k,v)->newMap.addToValue(k>>>shr,v));
        if(newMap.size()>=maxSize){
          SHR+=delta;
        }
      }
      SHR++;
      newMap.clear();
      final int shr = SHR;
      hashMap.forEachKeyValue((k,v)->newMap.addToValue(k>>>shr,v));
    }
    remainingBits -= SHR;
    hashMap = newMap;
  }

//  public void insert(long num, long freq) {
//    num >>>= bitsOfValue - remainingBits;
//    hashMap.addToValue(num,freq);
//    if (hashMap.size() == maxSize)
//      rebuild();
//  }
  private long doubleToLongBits(double data){
    long longBits = Double.doubleToLongBits(data) + (1L << 63);
    return data >= 0d ? longBits : longBits ^ 0x7FF0000000000000L;
  }

  public void insertDouble(double data) {
    long longBits = doubleToLongBits(data);
    longBits >>>= bitsOfValue - remainingBits;
    hashMap.addToValue(longBits,1);
    if (hashMap.size() == maxSize)
      rebuild();
  }
  public void insertLongBits(long longBits, long freq) {
    hashMap.addToValue(longBits,freq);
    if (hashMap.size() == maxSize)
      rebuild();
  }

  public int getRemainingBits() {
    return remainingBits;
  }

  public void merge(GSHashMapForStat b){
    if(b.remainingBits<remainingBits){
      final int delta = remainingBits - b.remainingBits;
      MutableLongLongMap newMap = new LongLongHashMap(maxSize);
      hashMap.forEachKeyValue((k,v)->newMap.addToValue(k>>>delta,v));
      hashMap = newMap;
      remainingBits-=delta;

      for(LongLongPair p:hashMap.keyValuesView()) {
        insertLongBits(p.getOne(), p.getTwo());
      }
    }else{
      final int delta = b.remainingBits - remainingBits;
      for(LongLongPair p:hashMap.keyValuesView()) {
        insertLongBits(p.getOne()>>>delta, p.getTwo());
      }
    }
  }

}
