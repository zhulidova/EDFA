function X = sigma_lum_param() 
% X(:,2) = B_e*10^(-21)
% X(:,3) = CrossEmTemp0*10^(-25)
X = [1445.8 -128.1793    0.0000
    1446.0 -105.1247    0.0000
    1446.3  -91.3784    0.0000
    1446.5  -83.1697    0.0000
    1446.7  -79.3944    0.0000
    1447.0  -74.8735    0.0000
    1447.2  -69.8495    0.0000
    1447.5  -65.4516    0.0000
    1447.8  -59.0505    0.0000
    1448.0  -53.7394    0.0000
    1448.3  -51.3741    0.0000
    1448.5  -48.9984    0.0000
    1448.7  -45.8056    0.0000
    1449.0  -43.1554    0.0000
    1449.2  -41.5587    0.0000
    1449.5  -40.3747    0.0000
    1449.8  -38.9073    0.0000
    1450.0  -37.1957    0.0000
    1450.3  -34.4350    0.0000
    1450.5  -31.6730    0.0000
    1450.7  -30.9597    0.0000
    1451.0  -29.9135    0.0000
    1451.2  -27.9720    0.0000
    1451.5  -26.3213    0.0001
    1451.8  -25.9497    0.0001
    1452.0  -25.4029    0.0001
    1452.3  -22.9237    0.0002
    1452.5  -20.9922    0.0003
    1452.7  -20.2912    0.0004
    1453.0  -19.8806    0.0005
    1453.2  -18.7214    0.0006
    1453.5  -17.7051    0.0008
    1453.8  -16.9592    0.0010
    1454.0  -16.1942    0.0013
    1454.3  -16.1019    0.0013
    1454.5  -15.8364    0.0014
    1454.7  -14.7659    0.0019
    1455.0  -13.8763    0.0025
    1455.2  -13.7004    0.0027
    1455.5  -13.6829    0.0027
    1455.8  -13.1010    0.0032
    1456.0  -12.7092    0.0036
    1456.3  -12.2123    0.0042
    1456.5  -11.6514    0.0050
    1456.7  -11.0507    0.0059
    1457.0  -10.4649    0.0069
    1457.2   -9.8205    0.0084
    1457.5   -9.1290    0.0103
    1457.8   -8.9413    0.0110
    1458.0   -8.6479    0.0120
    1458.3   -8.2315    0.0136
    1458.5   -7.7530    0.0157
    1458.7   -7.3222    0.0178
    1459.0   -6.6880    0.0212
    1459.2   -6.6752    0.0217
    1459.5   -6.6948    0.0219
    1459.8   -6.5467    0.0233
    1460.0   -6.3790    0.0248
    1460.3   -5.7254    0.0296
    1460.5   -5.3085    0.0333
    1460.7   -4.8815    0.0381
    1461.0   -4.4083    0.0440
    1461.2   -4.4644    0.0440
    1461.5   -4.5708    0.0435
    1461.8   -4.2107    0.0485
    1462.0   -4.0680    0.0511
    1462.3   -3.8279    0.0555
    1462.5   -3.6075    0.0599
    1462.7   -3.3972    0.0637
    1463.0   -3.1384    0.0686
    1463.2   -3.0084    0.0728
    1463.5   -2.6290    0.0822
    1463.8   -2.5279    0.0855
    1464.0   -2.5055    0.0873
    1464.3   -2.2741    0.0940
    1464.5   -2.0585    0.1005
    1464.7   -1.8863    0.1068
    1465.0   -1.5784    0.1174
    1465.2   -1.4600    0.1225
    1465.5   -1.5600    0.1212
    1465.8   -1.2769    0.1328
    1466.0   -0.9787    0.1461
    1466.3   -0.8663    0.1519
    1466.5   -0.8444    0.1541
    1466.7   -0.5335    0.1698
    1467.0   -0.2598    0.1852
    1467.2   -0.3011    0.1867
    1467.5   -0.1235    0.1989
    1467.8   -0.1952    0.1976
    1468.0   -0.2177    0.1985
    1468.3   -0.0640    0.2111
    1468.5   -0.0488    0.2166
    1468.8    0.0484    0.2231
    1469.0    0.1649    0.2303
    1469.2    0.3011    0.2430
    1469.5    0.5193    0.2620
    1469.8    0.5853    0.2706
    1470.0    0.7653    0.2869
    1470.3    0.7181    0.2856
    1470.5    0.6716    0.2836
    1470.8    1.0148    0.3175
    1471.0    1.2457    0.3451
    1471.2    1.1825    0.3432
    1471.5    1.3339    0.3596
    1471.8    1.5225    0.3821
    1472.0    1.6920    0.4042
    1472.3    1.7043    0.4107
    1472.5    1.6002    0.4051
    1472.8    1.8128    0.4335
    1473.0    1.9414    0.4538
    1473.2    2.1289    0.4837
    1473.5    2.2219    0.5044
    1473.8    2.1520    0.4993
    1474.0    2.2736    0.5181
    1474.3    2.3338    0.5352
    1474.5    2.4149    0.5573
    1474.8    2.5106    0.5755
    1475.0    2.6375    0.5988
    1475.2    2.7366    0.6200
    1475.5    2.7503    0.6275
    1475.8    2.7636    0.6415
    1476.0    2.6087    0.6293
    1476.3    2.8946    0.6782
    1476.5    2.9786    0.6957
    1476.8    3.0708    0.7255
    1477.0    3.0652    0.7373
    1477.2    3.0811    0.7457
    1477.5    3.1876    0.7702
    1477.8    3.1059    0.7626
    1478.0    3.0390    0.7581
    1478.3    3.1680    0.7987
    1478.5    3.2838    0.8389
    1478.8    3.2345    0.8294
    1479.0    3.2860    0.8394
    1479.2    3.2420    0.8472
    1479.5    3.3092    0.8775
    1479.8    3.2860    0.8753
    1480.0    3.2526    0.8709
    1480.3    3.3385    0.8994
    1480.5    3.3454    0.9120
    1480.8    3.3990    0.9400
    1481.0    3.4508    0.9685
    1481.2    3.3915    0.9558
    1481.5    3.4332    0.9681
    1481.8    3.5431    1.0128
    1482.0    3.4735    1.0138
    1482.3    3.5436    1.0360
    1482.5    3.6138    1.0587
    1482.8    3.6561    1.0796
    1483.0    3.6775    1.0952
    1483.2    3.6887    1.1146
    1483.5    3.7058    1.1328
    1483.8    3.7897    1.1567
    1484.0    3.7945    1.1570
    1484.3    3.7947    1.1767
    1484.5    3.6901    1.1677
    1484.8    3.7312    1.1840
    1485.0    3.6850    1.1764
    1485.2    3.8621    1.2398
    1485.5    3.8631    1.2512
    1485.8    3.9228    1.2863
    1486.0    3.9857    1.3236
    1486.3    3.8850    1.2901
    1486.5    3.8827    1.2901
    1486.8    3.8966    1.3156
    1487.0    3.9197    1.3448
    1487.2    4.0730    1.3981
    1487.5    4.0517    1.3924
    1487.8    4.0656    1.4057
    1488.0    4.0654    1.4177
    1488.3    4.0995    1.4490
    1488.5    4.1329    1.4807
    1488.8    4.1214    1.4750
    1489.0    4.1890    1.4999
    1489.2    4.1909    1.5243
    1489.5    4.1753    1.5385
    1489.8    4.1774    1.5492
    1490.0    4.1775    1.5553
    1490.3    4.1882    1.5684
    1490.5    4.2074    1.5852
    1490.8    4.1925    1.5965
    1491.0    4.2632    1.6400
    1491.2    4.2009    1.6150
    1491.5    4.1302    1.5832
    1491.8    4.1364    1.6112
    1492.0    4.2359    1.6758
    1492.3    4.2639    1.6940
    1492.5    4.2776    1.7063
    1492.8    4.2795    1.7128
    1493.0    4.2773    1.7221
    1493.2    4.2985    1.7540
    1493.5    4.3211    1.7917
    1493.8    4.2344    1.7531
    1494.0    4.2249    1.7459
    1494.3    4.2093    1.7616
    1494.5    4.2787    1.8128
    1494.8    4.2768    1.8227
    1495.0    4.2749    1.8327
    1495.2    4.2599    1.8324
    1495.5    4.1716    1.7973
    1495.8    4.2992    1.8817
    1496.0    4.2620    1.8920
    1496.3    4.2530    1.8845
    1496.5    4.2495    1.8794
    1496.8    4.2582    1.9042
    1497.0    4.1699    1.8861
    1497.2    4.1734    1.9031
    1497.5    4.0854    1.8753
    1497.8    4.1049    1.8861
    1498.0    4.0452    1.8583
    1498.3    4.0515    1.8912
    1498.5    4.0579    1.9246
    1498.8    3.9737    1.8890
    1499.0    3.9737    1.8906
    1499.2    3.8759    1.8587
    1499.5    3.8593    1.8667
    1499.8    3.8479    1.8813
    1500.0    3.9274    1.9412
    1500.3    3.8616    1.9098
    1500.5    3.7957    1.8788
    1500.8    3.7330    1.8762
    1501.0    3.7621    1.9189
    1501.2    3.7541    1.9173
    1501.5    3.7366    1.9159
    1501.8    3.7699    1.9434
    1502.0    3.7917    1.9657
    1502.3    3.7123    1.9515
    1502.5    3.7077    1.9757
    1502.8    3.6226    1.9285
    1503.0    3.6118    1.9195
    1503.2    3.5420    1.9143
    1503.5    3.6466    1.9929
    1503.8    3.5378    1.9488
    1504.0    3.5142    1.9484
    1504.3    3.4484    1.9205
    1504.5    3.3797    1.8965
    1504.8    3.4447    1.9553
    1505.0    3.4271    1.9778
    1505.2    3.3480    1.9391
    1505.5    3.3581    1.9455
    1505.7    3.2853    1.9337
    1506.0    3.2931    1.9629
    1506.3    3.2849    1.9744
    1506.5    3.2808    1.9832
    1506.8    3.2702    1.9849
    1507.0    3.1733    1.9422
    1507.2    3.2023    1.9867
    1507.5    3.2251    2.0289
    1507.7    3.1460    1.9938
    1508.0    3.1515    1.9983
    1508.3    3.1466    2.0168
    1508.5    3.1380    2.0338
    1508.8    3.0544    2.0156
    1509.0    3.0514    2.0402
    1509.2    3.0621    2.0475
    1509.5    2.8964    1.9676
    1509.7    3.0233    2.0646
    1510.0    3.0595    2.1158
    1510.3    3.0512    2.1239
    1510.5    2.9604    2.0917
    1510.8    2.9518    2.1034
    1511.0    2.9382    2.1125
    1511.2    2.9518    2.1526
    1511.5    2.9602    2.1904
    1511.7    2.8977    2.1614
    1512.0    2.9983    2.2202
    1512.3    2.9870    2.2477
    1512.5    2.9759    2.2756
    1512.8    2.9009    2.2548
    1513.0    2.9083    2.2771
    1513.2    2.9337    2.3042
    1513.5    2.9601    2.3326
    1513.7    2.9482    2.3726
    1514.0    3.0203    2.4610
    1514.3    2.9591    2.4353
    1514.5    2.8977    2.4042
    1514.8    2.8971    2.4352
    1515.0    2.8976    2.4733
    1515.2    2.9211    2.5256
    1515.5    2.8640    2.5253
    1515.7    2.7933    2.4926
    1516.0    2.8015    2.5058
    1516.3    2.7963    2.5535
    1516.5    2.7885    2.6006
    1516.8    2.7145    2.5778
    1517.0    2.8083    2.6627
    1517.2    2.7341    2.6514
    1517.5    2.6627    2.6357
    1517.7    2.7530    2.7527
    1518.0    2.7576    2.8115
    1518.3    2.5860    2.7105
    1518.5    2.5008    2.6726
    1518.8    2.6736    2.8518
    1519.0    2.6784    2.9195
    1519.2    2.6031    2.9134
    1519.5    2.6069    2.9610
    1519.7    2.5331    2.9419
    1520.0    2.4608    2.9170
    1520.3    2.4788    3.0104
    1520.5    2.5007    3.1023
    1520.8    2.4134    3.0791
    1521.0    2.4076    3.1143
    1521.2    2.3349    3.1168
    1521.5    2.3461    3.1880
    1521.7    2.2795    3.2187
    1522.0    2.2144    3.2433
    1522.3    2.2131    3.2836
    1522.5    2.2904    3.3851
    1522.8    2.2954    3.4815
    1523.0    2.2984    3.5789
    1523.2    2.1427    3.5054
    1523.5    2.0685    3.5071
    1523.7    2.0746    3.5731
    1524.0    2.1609    3.7081
    1524.3    2.1601    3.7466
    1524.5    1.9218    3.5918
    1524.8    1.8459    3.5660
    1525.0    1.6089    3.3947
    1525.2    1.1311    3.0570
    1525.5    0.6566    2.7485
    1525.7   -0.1573    2.2894
    1526.0   -0.2366    2.2657
    1526.3   -0.1619    2.3178
    1526.5   -0.0067    2.4212
    1526.8    0.4827    2.7657
    1527.0    0.3999    2.7219
    1527.2    0.4027    2.7325
    1527.5    0.4888    2.7966
    1527.7    0.7309    2.9743
    1528.0    0.5726    2.8766
    1528.3    0.6446    2.9476
    1528.5    0.6321    2.9479
    1528.8    0.0757    2.6003
    1529.0   -0.0769    2.5066
    1529.2   -0.2358    2.4297
    1529.5   -0.0796    2.5268
    1529.7   -0.0797    2.5227
    1530.0   -0.0735    2.5347
    1530.3    0.1692    2.6895
    1530.5    0.1693    2.6852
    1530.8    0.2470    2.7425
    1531.0    0.2443    2.7364
    1531.3    0.0806    2.6307
    1531.5    0.2415    2.7265
    1531.7    0.3144    2.7552
    1532.0    0.2297    2.6785
    1532.3   -0.3459    2.2967
    1532.5   -0.9199    1.9747
    1532.8   -0.6786    2.0703
    1533.0   -0.9208    1.9135
    1533.3   -1.0077    1.8504
    1533.5   -1.1754    1.7563
    1533.7   -1.0903    1.7763
    1534.0   -1.0881    1.7622
    1534.3   -0.4466    2.0359
    1534.5   -0.1296    2.1665
    1534.8   -0.3758    2.0144
    1535.0   -0.2935    2.0404
    1535.3   -0.2202    2.0572
    1535.5   -0.3155    1.9525
    1535.7   -0.5676    1.7998
    1536.0   -0.9785    1.5919
    1536.3   -0.7351    1.6815
    1536.5   -0.8207    1.6341
    1536.8   -1.1425    1.5074
    1537.0   -1.1464    1.4999
    1537.3   -1.3017    1.4449
    1537.5   -1.2387    1.4597
    1537.7   -0.8247    1.6157
    1538.0   -0.5735    1.7182
    1538.3   -0.2469    1.8459
    1538.5   -0.2459    1.8434
    1538.8   -0.2534    1.8414
    1539.0   -0.3484    1.7944
    1539.3   -0.4238    1.7608
    1539.5   -0.3466    1.7854
    1539.7   -0.2694    1.8146
    1540.0   -0.3687    1.7613
    1540.3   -0.2096    1.8244
    1540.5    0.0363    1.9330
    1540.8    0.0285    1.9260
    1541.0    0.0255    1.9213
    1541.3    0.2646    2.0366
    1541.5    0.1024    1.9579
    1541.7    0.2675    2.0407
    1542.0    0.2745    2.0556
    1542.3    0.1864    2.0060
    1542.5    0.1034    1.9600
    1542.8    0.1064    1.9587
    1543.0    0.1944    2.0153
    1543.3    0.1034    1.9652
    1543.5   -0.0658    1.8774
    1543.7   -0.0659    1.8656
    1544.0    0.1804    1.9814
    1544.3   -0.0740    1.8605
    1544.5   -0.0780    1.8646
    1544.8   -0.0019    1.8947
    1545.0    0.1601    1.9593
    1545.3   -0.0024    1.8749
    1545.5   -0.1566    1.7933
    1545.7   -0.1555    1.7995
    1546.0   -0.1535    1.7977
    1546.3   -0.0819    1.8196
    1546.5   -0.0115    1.8327
    1546.8   -0.1002    1.7926
    1547.0    0.1435    1.9028
    1547.3   -0.0928    1.7939
    1547.5   -0.0818    1.7910
    1547.7   -0.0205    1.8041
    1548.0   -0.2036    1.7174
    1548.3    0.1490    1.8750
    1548.5   -0.0759    1.7737
    1548.8    0.0813    1.8361
    1549.0    0.0702    1.8150
    1549.3    0.0874    1.8159
    1549.5    0.1860    1.8513
    1549.7    0.3442    1.9263
    1550.0    0.5814    2.0413
    1550.3    0.2516    1.8791
    1550.5    0.2429    1.8585
    1550.8    0.1601    1.8165
    1551.0    0.3215    1.8920
    1551.3    0.3351    1.9091
    1551.5    0.2664    1.8764
    1551.7    0.5017    1.9825
    1552.0   -0.4876    1.5480
    1552.3    0.2405    1.8381
    1552.5    0.0705    1.7646
    1552.8    0.3125    1.8763
    1553.0    0.3859    1.9050
    1553.3    0.3046    1.8578
    1553.5    0.3850    1.8853
    1553.7    0.3155    1.8567
    1554.0    0.0821    1.7563
    1554.3    0.3984    1.8999
    1554.5    0.4778    1.9229
    1554.8    0.4001    1.8813
    1555.0    0.4086    1.8824
    1555.3    0.4093    1.8935
    1555.5    0.6570    2.0167
    1555.7    0.5678    1.9625
    1556.0    0.6543    1.9992
    1556.3    0.4051    1.8709
    1556.5    0.6713    2.0024
    1556.8    0.5131    1.9276
    1557.0    0.6822    2.0063
    1557.3    0.6803    1.9928
    1557.5    0.9276    2.1021
    1557.7    0.5927    1.9334
    1558.0    0.8363    2.0574
    1558.3    0.7469    2.0070
    1558.5    0.7459    1.9987
    1558.8    0.6504    1.9334
    1559.0    0.4037    1.8022
    1559.3    0.7260    1.9569
    1559.5    0.9667    2.0805
    1559.7    0.6314    1.8995
    1560.0    0.6267    1.8768
    1560.3    0.7003    1.9008
    1560.5    0.4440    1.7708
    1560.8    0.6020    1.8428
    1561.0    0.6784    1.8819
    1561.3    0.9949    2.0121
    1561.5    1.3978    2.2055
    1561.7    1.6451    2.3319
    1562.0    1.4042    2.1968
    1562.3    1.5675    2.2727
    1562.5    1.7370    2.3549
    1562.8    1.6344    2.2897
    1563.0    2.1147    2.5598
    1563.3    2.1214    2.5597
    1563.5    1.1401    1.9831
    1563.7    0.7838    1.7788
    1564.0    1.0321    1.8721
    1564.3    1.2696    1.9699
    1564.5    1.1180    1.8998
    1564.8    0.7733    1.7465
    1565.0    0.6752    1.7000
    1565.3    1.2493    1.9388
    1565.5    1.5701    2.0652
    1565.7    1.5932    2.0632
    1566.0    1.6932    2.1029
    1566.3    1.6010    2.0398
    1566.5    1.6726    2.0597
    1566.8    1.5804    1.9745
    1567.0    1.8174    2.0623
    1567.3    1.9003    2.0936
    1567.5    1.9832    2.1304
    1567.7    1.8900    2.0517
    1568.0    2.2984    2.2516
    1568.2    2.3729    2.2765
    1568.5    1.6967    1.9113
    1568.8    1.7032    1.9112
    1569.0    1.7037    1.9037
    1569.3    1.7712    1.9027
    1569.5    2.1724    2.0742
    1569.7    2.5908    2.2807
    1570.0    3.0149    2.5114
    1570.2    3.1071    2.5565
    1570.5    3.0315    2.4912
    1570.8    3.3367    2.6477
    1571.0    3.2277    2.5323
    1571.3    3.2112    2.5000
    1571.5    3.1141    2.4166
    1571.7    3.1273    2.4023
    1572.0    3.1333    2.3895
    1572.2    3.0251    2.2820
    1572.5    2.9094    2.1698
    1572.8    2.9063    2.1544
    1573.0    2.8218    2.0994
    1573.3    2.9091    2.1177
    1573.5    2.9927    2.1393
    1573.7    2.7968    2.0004
    1574.0    2.6827    1.9063
    1574.2    2.5977    1.8661
    1574.5    2.5960    1.8622
    1574.8    2.6080    1.8380
    1575.0    2.6160    1.8124
    1575.3    2.5110    1.7390
    1575.5    2.4861    1.7038
    1575.7    2.3779    1.6573
    1576.0    2.2691    1.6118
    1576.2    2.2373    1.5694
    1576.5    2.2769    1.5531
    1576.8    2.1798    1.5012
    1577.0    2.0777    1.4526
    1577.3    2.0642    1.4408
    1577.5    2.0384    1.4247
    1577.7    1.9355    1.3646
    1578.0    1.9185    1.3335
    1578.2    1.8549    1.3054
    1578.5    1.7814    1.2779
    1578.8    1.7512    1.2625
    1579.0    1.6336    1.2192
    1579.3    1.7290    1.2292
    1579.5    1.7381    1.2116
    1579.7    1.5449    1.1534
    1580.0    1.5287    1.1468
    1580.2    1.4280    1.1100
    1580.5    1.4061    1.0966
    1580.8    1.3132    1.0576
    1581.0    1.2964    1.0378
    1581.3    1.1772    1.0098
    1581.5    1.1516    1.0044
    1581.7    1.0517    0.9720
    1582.0    0.9514    0.9405
    1582.2    0.9321    0.9244
    1582.5    0.8366    0.8929
    1582.8    0.8209    0.8917
    1583.0    0.6314    0.8532
    1583.3    0.6331    0.8432
    1583.5    0.6191    0.8302
    1583.7    0.5179    0.8058
    1584.0    0.4957    0.7965
    1584.2    0.5040    0.7999
    1584.5    0.5192    0.8045
    1584.8    0.4358    0.7778
    1585.0    0.3694    0.7551
    1585.3    0.3221    0.7449
    1585.5    0.4429    0.7659
    1585.7    0.3114    0.7418
    1586.0    0.2593    0.7335
    1586.2    0.1755    0.7096
    1586.5    0.0891    0.6860
    1586.8   -0.0019    0.6700
    1587.0    0.0111    0.6706
    1587.3    0.0688    0.6810
    1587.5    0.0326    0.6748
    1587.7    0.0030    0.6617
    1588.0    0.0592    0.6620
    1588.2    0.0441    0.6613
    1588.5   -0.0501    0.6471
    1588.8    0.0243    0.6571
    1589.0    0.0904    0.6660
    1589.3   -0.0770    0.6337
    1589.5   -0.2691    0.5977
    1589.7   -0.2122    0.6085
    1590.0   -0.1599    0.6172
    1590.2   -0.1599    0.6148
    1590.5   -0.1487    0.6143
    1590.8   -0.2730    0.5893
    1591.0   -0.3942    0.5673
    1591.3   -0.3524    0.5751
    1591.5   -0.3034    0.5841
    1591.7   -0.3717    0.5703
    1592.0   -0.3504    0.5699
    1592.2   -0.5563    0.5365
    1592.5   -0.7707    0.5052
    1592.8   -0.6960    0.5154
    1593.0   -0.6125    0.5285
    1593.3   -0.6423    0.5208
    1593.5   -0.6726    0.5130
    1593.8   -0.6739    0.5094
    1594.0   -0.7418    0.4982
    1594.2   -0.6956    0.5059
    1594.5   -0.8201    0.4925
    1594.8   -0.8233    0.4862
    1595.0   -0.9263    0.4688
    1595.3   -0.8053    0.4817
    1595.5   -0.7553    0.4859
    1595.8   -0.9170    0.4686
    1596.0   -0.9078    0.4714
    1596.2   -0.9052    0.4661
    1596.5   -0.9944    0.4500
    1596.8   -0.9719    0.4518
    1597.0   -1.1085    0.4363
    1597.3   -1.0030    0.4488
    1597.5   -0.9850    0.4524
    1597.8   -0.9911    0.4454
    1598.0   -1.1064    0.4274
    1598.2   -1.1382    0.4241
    1598.5   -1.0776    0.4311
    1598.8   -1.1028    0.4286
    1599.0   -1.0580    0.4341
    1599.3   -0.9867    0.4354
    1599.5   -1.0959    0.4177
    1599.8   -0.9981    0.4289
    1600.0   -0.9142    0.4388
    1600.2   -0.8650    0.4434
    1600.5   -0.7182    0.4582
    1600.8   -0.8009    0.4433
    1601.0   -0.9649    0.4199
    1601.3   -0.8105    0.4386
    1601.5   -0.8522    0.4364
    1601.8   -0.7641    0.4433
    1602.0   -0.7595    0.4406
    1602.2   -0.8431    0.4270
    1602.5   -0.8430    0.4230
    1602.8   -0.7774    0.4329
    1603.0   -0.6065    0.4553
    1603.3   -0.6710    0.4441
    1603.5   -0.7221    0.4346
    1603.8   -0.7735    0.4264
    1604.0   -0.8346    0.4173
    1604.2   -0.8293    0.4201
    1604.5   -0.6534    0.4412
    1604.8   -0.7966    0.4212
    1605.0   -0.9129    0.4059
    1605.3   -0.9744    0.3966
    1605.5   -1.0229    0.3888
    1605.8   -1.1576    0.3787
    1606.0   -1.1910    0.3787
    1606.2   -1.1474    0.3779
    1606.5   -1.2143    0.3673
    1606.8   -1.2632    0.3612
    1607.0   -1.4257    0.3458
    1607.3   -1.3696    0.3523
    1607.5   -1.4218    0.3490
    1607.8   -1.4745    0.3402
    1608.0   -1.6478    0.3215
    1608.2   -1.6348    0.3218
    1608.5   -1.6060    0.3234
    1608.8   -1.7294    0.3144
    1609.0   -1.7424    0.3146
    1609.3   -1.6212    0.3202
    1609.5   -1.6072    0.3170
    1609.8   -1.5636    0.3203
    1610.0   -1.6048    0.3165
    1610.2   -1.6119    0.3158
    1610.5   -1.6292    0.3143
    1610.8   -1.8209    0.2962
    1611.0   -1.8445    0.2909
    1611.3   -1.9019    0.2872
    1611.5   -2.0580    0.2771
    1611.8   -1.9809    0.2824
    1612.0   -1.9133    0.2872
    1612.2   -2.2428    0.2607
    1612.5   -2.4644    0.2434
    1612.8   -2.3273    0.2525
    1613.0   -2.2077    0.2608
    1613.3   -2.2378    0.2588
    1613.5   -2.3645    0.2504
    1613.8   -2.4670    0.2408
    1614.0   -2.4690    0.2369
    1614.2   -2.5406    0.2336
    1614.5   -2.5218    0.2353
    1614.8   -2.5975    0.2307
    1615.0   -2.6552    0.2273
    1615.3   -2.5615    0.2293
    1615.5   -2.6542    0.2210
    1615.8   -2.6177    0.2243
    1616.0   -2.6605    0.2230
    1616.2   -2.7517    0.2169
    1616.5   -2.9157    0.2076
    1616.8   -3.0135    0.2000
    1617.0   -3.0944    0.1936
    1617.3   -3.1509    0.1916
    1617.5   -3.1880    0.1907
    1617.8   -3.2390    0.1872
    1618.0   -3.2711    0.1847
    1618.2   -3.3754    0.1780
    1618.5   -3.4618    0.1724
    1618.8   -3.4879    0.1716
    1619.0   -3.6053    0.1668
    1619.3   -3.5571    0.1681
    1619.5   -3.6328    0.1644
    1619.8   -3.7316    0.1587
    1620.0   -3.7631    0.1554
    1620.2   -3.8297    0.1532
    1620.5   -3.9088    0.1505
    1620.8   -3.8770    0.1508
    1621.0   -3.7302    0.1553
    1621.3   -3.7577    0.1524
    1621.5   -3.6951    0.1531
    1621.8   -3.8685    0.1469
    1622.0   -3.9499    0.1444
    1622.2   -3.9738    0.1426
    1622.5   -4.0073    0.1402
    1622.8   -4.0776    0.1363
    1623.0   -4.1617    0.1320
    1623.3   -4.2329    0.1300
    1623.5   -4.4911    0.1223
    1623.8   -4.4322    0.1233
    1624.0   -4.5155    0.1199
    1624.2   -4.5468    0.1175
    1624.5   -4.4849    0.1179
    1624.8   -4.6033    0.1150
    1625.0   -4.8261    0.1092
    1625.3   -4.9175    0.1054
    1625.5   -4.9966    0.1021
    1625.8   -4.9668    0.1020
    1626.0   -4.8578    0.1039
    1626.2   -4.8502    0.1045
    1626.5   -4.7215    0.1084
    1626.8   -5.0115    0.0999
    1627.0   -5.3033    0.0917
    1627.3   -5.1857    0.0935
    1627.5   -5.2967    0.0901
    1627.8   -5.3285    0.0897
    1628.0   -5.2373    0.0920
    1628.2   -5.2215    0.0913
    1628.5   -5.2043    0.0905
    1628.8   -5.2616    0.0885
    1629.0   -5.4940    0.0828
    1629.3   -5.3120    0.0868
    1629.5   -5.1024    0.0917
    1629.8   -5.4276    0.0839
    1630.0   -5.6613    0.0783
    1630.2   -5.5101    0.0804
    1630.5   -5.5626    0.0784
    1630.8   -5.5577    0.0785
    1631.0   -5.6493    0.0768
    1631.3   -5.7574    0.0743
    1631.5   -5.8796    0.0713
    1631.8   -6.0474    0.0677
    1632.0   -6.2306    0.0639
    1632.2   -5.9467    0.0688
    1632.5   -5.8053    0.0714
    1632.7   -6.1728    0.0646
    1633.0   -6.4779    0.0592
    1633.3   -6.6764    0.0556
    1633.5   -6.8943    0.0517
    1633.8   -6.6251    0.0558
    1634.0   -6.3912    0.0596
    1634.2   -6.3704    0.0591
    1634.5   -6.4904    0.0564
    1634.7   -6.6030    0.0547
    1635.0   -6.8266    0.0514
    1635.3   -6.7727    0.0519
    1635.5   -6.9305    0.0497
    1635.8   -6.7179    0.0519
    1636.0   -6.6711    0.0519
    1636.2   -6.9268    0.0482
    1636.5   -7.2674    0.0438
    1636.7   -7.0449    0.0465
    1637.0   -6.9581    0.0477
    1637.3   -7.1657    0.0447
    1637.5   -7.5141    0.0404
    1637.8   -7.5141    0.0401
    1638.0   -7.5141    0.0397
    1638.2   -7.4656    0.0400
    1638.5   -7.3181    0.0413
    1638.7   -7.5088    0.0390
    1639.0   -7.6803    0.0370
    1639.3   -7.8168    0.0354
    1639.5   -7.8321    0.0348
    1639.8   -7.9145    0.0342
    1640.0   -7.9171    0.0344
    1640.2   -8.1113    0.0323
    1640.5   -8.4284    0.0294
    1640.7   -8.3011    0.0300
    1641.0   -8.3575    0.0292
    1641.3   -8.2785    0.0299
    1641.5   -8.2587    0.0301
    1641.8   -8.5657    0.0277
    1642.0   -8.9960    0.0246
    1642.2   -9.0449    0.0240
    1642.5   -9.0521    0.0236
    1642.7   -9.1530    0.0231
    1643.0   -9.1749    0.0229
    1643.3   -9.3880    0.0215
    1643.5   -9.5807    0.0203
    1643.8   -9.5247    0.0204
    1644.0   -9.7207    0.0192
    1644.2   -9.8234    0.0187
    1644.5   -9.9738    0.0179
    1644.7  -10.3246    0.0163
    1645.0  -10.4962    0.0155
    1645.3  -10.5308    0.0151
    1645.5  -10.6907    0.0143
    1645.8  -10.8185    0.0138
    1646.0  -10.9215    0.0135
    1646.2  -11.2616    0.0123
    1646.5  -11.6068    0.0112
    1646.7  -11.8586    0.0103
    1647.0  -12.0626    0.0096
    1647.3  -12.0095    0.0098
    1647.5  -11.8503    0.0102
    1647.8  -12.1226    0.0095];
end 