# Exact-Quantile

Branch in Apache IoTDB: https://github.com/apache/iotdb/tree/research/exact-quantile

Three open datasets from Kaggle used in experiments . 


	Bitcoin Dataset:
		https://www.kaggle.com/datasets/shiheyingzhe/bitcoin-transaction-data-from-2009-to-2018
		https://userscloud.com/
		
	Thruster Dataset:
		https://www.kaggle.com/datasets/patrickfleith/spacecraft-thruster-firing-tests-dataset
		https://userscloud.com/
		
	Taxi Dataset:
		https://www.kaggle.com/datasets/sunnets/taxipredictiondata8m
		https://userscloud.com/





### Out-of-database test in paper:

The test code is in src/main/java.

The verification experiments: ***Fig. 12, 13, 14*** and the "Theoretical" results in ***Fig 15***.

​		Please download the 3 dataset and store in the root directory of this project.

**Fig 12, 15: **

```
	MainForExactLazyPrioriKLLShowFixedPrF
```

**Fig 13:**

```
	MainForVerifyFilterSize.java
```

**Fig 14:**

```
	MainForExactLazyPrioriKLLShowDistributionOfFilterSize.java.
```











### In-database test:

The observed time cost and passes are recorded after executing queries in IoTDB server.

​		**Part 1	Get IoTDB with exact-quantile:**

​			Clone https://github.com/apache/iotdb/tree/research/exact-quantile , then compile IoTDB server and client by 

				mvn spotless:apply
				mvn clean package -pl cli -am -Dmaven.test.skip=trueD:
​		**Part 2	Run IoTDB Server :**

​			We need to check the configuration file iotdb\server\target\iotdb-server-0.13.0-SNAPSHOT\conf\iotdb-engine.properties 

​			The main parameters needed to be confirmed are as follows:

| name in configuration file | meaning                                                |
| -------------------------- | ------------------------------------------------------ |
| synopsis_size_in_byte      | Size of the pre-computed sketch size for each Page     |
| enable_synopsis            | Whether to collect pre-computed sketches in ingestion. |

​			Then we can run server by 

				.\iotdb\server\target\iotdb-server-0.13.0-SNAPSHOT\sbin\start-server.bat
​				or

```
			.\iotdb\server\target\iotdb-server-0.13.0-SNAPSHOT\sbin\start-server.sh
```

​			It is required to <u>**restart**</u> server for any new configuration to take effect.



​		**Part 3	Ingest Data Set**

​				Download the three open datasets from the link above, and ingest them by running SessionSyn.py:

​				To ingest data without updates, it's recommended to execute iotdb\session\src\test\java\org\apache\iotdb\session\InsertCsvDataIT.java

​				To ingest data with updates for ***Fig. 21***, it's recommended to execute iotdb\session\src\test\java\org\apache\iotdb\session\ExactQuantile_InsertUpdateData.java

​				Please look into the java code to determine the path of dataset and its schema in IoTDB.



​		**Part 4 TEST**

​				Methods are different aggregation functions deployed in IoTDB.

​				One can check all aggregate functions or all ingested timeseries (datasets) by executing "show functions" or "show timeseries" in IoTDB Client.



​				**Test in IoTDB client**

​					One can run client and execute any query by 

				.\iotdb\cli\target\iotdb-cli-0.13.0-SNAPSHOT\sbin\start-cli.bat -h 127.0.0.1 -p 6667 -u root -pw root

​				or

```
			.\iotdb\cli\target\iotdb-cli-0.13.0-SNAPSHOT\sbin\start-cli.sh -h 127.0.0.1 -p 6667 -u root -pw root
```

​						Note that the client will show the time cost of executed query statement.



SQL statements for querying single quantile (***Fig. 15, 16, 17, 20, 21***):

```
Det-MRL: exact_quantile_mrl(s0,'quantile'='0.5','memory'='1024KB','return_type'='iteration_num','merge_buffer_ratio'='0') from ... where time ......

Det-t-digest:
exact_quantile_tdigest(s0,'quantile'='0.5','memory'='1024KB','return_type'='iteration_num','param'='1') from ... where time ......

Det-DDSketch:
exact_quantile_ddsketch_positive(s0,'quantile'='0.5','memory'='1024KB','return_type'='iteration_num') from ... where time ......

Rand-KLL:
exact_quantile_pr_kll_post_best_pr(s0,'quantile'='0.5','memory'='1024KB','return_type'='iteration_num','merge_buffer_ratio'='0') from ... where time ......

Rand-Page-KLL:
exact_quantile_pr_kll_post_best_pr(s0,'multi_quantiles'='32','memory'='1024KB','return_type'='iteration_num','merge_buffer_ratio'='5') from DataWithPreComputedSynopses from ... where time ......

```

*'return_type'='iteration_num'* means that the return value of the statement is the number of passes.

​	It can be changed to like *'return_type'='value'*  to show the computed quantile.

The amount of data (**N**) is controlled by the where clause like " where time>=0 and time<100 ".

The memory budget (**Memory**) is controlled by the parameter 'memory'='...'

The queried quantile is controlled by the parameter 'quantile'='...'

​			

SQL statements for querying multiple quantiles (***Fig. 18***):

```
Det-MRL: exact_multi_quantiles_mrl(s0,'multi_quantiles'='32','memory'='1024KB','return_type'='iteration_num','merge_buffer_ratio'='0') from ... where time ......

Det-t-digest:
exact_multi_quantiles_tdigest(s0,'multi_quantiles'='32','memory'='1024KB','return_type'='iteration_num','param'='1') from ... where time ......

Det-DDSketch:
exact_multi_quantiles_ddsketch_positive(s0,'multi_quantiles'='32','memory'='1024KB','return_type'='iteration_num') from ... where time ......

Rand-KLL:
exact_multi_quantiles_pr_kll_post_best_pr(s0,'multi_quantiles'='32','memory'='1024KB','return_type'='iteration_num','merge_buffer_ratio'='0') from ... where time ......

Rand-Page-KLL:
exact_multi_quantiles_pr_kll_post_best_pr(s0,'multi_quantiles'='32','memory'='1024KB','return_type'='iteration_num','merge_buffer_ratio'='5')  from ...(DataWithPreComputedSynopses) where time ...... 

```

The number of queried quantiles (**Quantile quantity**) is controlled by the parameter 'multi_quantiles'='...'





Ingestion time record for ***Fig. 19(b)***: 

​	iotdb\session\src\test\java\org\apache\iotdb\session\InsertCsvDataIT.java



