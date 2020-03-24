# Вычисление коэффициентов светосбора для программы моделирования счетчиков АЧС

Для вычисления коэффициентов светосбора данные полученные на космических мюонах разбиваются на временные промежутки 1.5 месяца. <br />
light_collection - Код для получения коэффициентов светосбора для каждой из области разбиения счетчика (28,58,68 - областей) - усреднение по каждому из пяти типов счетчиков. <br /> 
Для запуска (определение коэффициентов для 1 типа счетчиков для 10 временного диапазона):
```
/light_collection 10 1 13122018 >> /dev/null &   
```
Для запуска на batch: <br />
```
qsub -shell n -b y -V -N reg1 -cwd ./light_collection 1 1 11022019
```
light_col_coef_cnt - Вычисление среднего поправочного коэффициента, который в моделировнии домножается на коэффициент светосбора, для каждого счетчика по всем разбиением. <br />
Для запуска:
``` 
/light_col_coef_cnt 16 1 12122018
```

# Временной интервал для получения коэффициентов неоднородности из заходов от космических мюонов:

1) 1 область <br />
   29.04.2014 - поднято поле в КЕДРе. <br />
   2914085 событий <br />     
   старт - RUN=19635    дата: 23.05.2014 - поле 6 кГс       ->   cosm_runs_jun14_1.root  - 2106 - номер калибровки- 22/05/2014 <br />
   cosm_runs_jun14_2.root  <br />
   cosm_runs_jun14_3.root  <br />
   cosm_runs_jun14_4.root  <br />
   cosm_runs_jun14_5.root  <br />
   cosm_runs_jun14_6.root  <br />
   cosm_runs_jun14_7.root  <br />
   RUN=19803    дата: 11.06.2014                   ->   cosm_runs_jun14_8.root  <br />
   cosm_runs_jun14_9.root  <br />                                  
   RUN=19810    дата: 12.06.2014                   ->   cosm_runs_jun14_9_1.root <br />
   3738414 событий в диапазоне заходов    <br />      
   финиш - RUN=19823    дата: 13.06.2014              ->   cosm_runs_jun14_10.root <br /> 
   cosm_runs_jun14_11.root <br />
   cosm_runs_jun14_12.root <br />
   cosm_runs_jun14_13.root <br />   
   cosm_runs_jun14_14.root <br />
   финиш - RUN=19854    дата: 16.06.2014 - поле 6 кГс  ->   cosm_runs_jun14_15.root <br />          
   1668976 событий в захоах:<br />
   старт- RUN=19862    дата: 18.06.2014 -     ->   cosm_runs_jun14_16.root - 2197 - номер калибровки - 16/06/14 <br />
   старт - RUN=19911    дата: 22.06.2014      ->   cosm_runs_jun14_17.root <br />
   финиш - RUN=19923    дата: 23.06.2014            <br />
   старт - RUN=19925    дата: 23.06.2014      ->   cosm_runs_jun14_18.root <br />              
   финиш - RUN=19929    дата: 24.06.2014              <br />
   старт - RUN=19930    дата: 24.06.2014      ->   cosm_runs_jun14_19.root <br />
   финиш - RUN=19936    дата: 25.06.2014 <br />
                  
   Заходы не включены в 1 область:        <br />            
   26.06.2014 - Поле 7000 Гс <br /> 
   27.06.2014 - Поле 7200 Гс <br />
   1.07.2014 - Плановый сброс магнитного поля. <br />
                                           
   14.09.2014 - Поле 2.8 кГс. <br />

   старт - RUN=19988    дата: 15.09.2014      <br />
   финиш - RUN=19995    дата: 15.09.2014      <br />
   старт - RUN=19996    дата: 15.09.2014      <br />
   финиш - RUN=19999    дата: 15.09.2014     <br />
                           
   18.09.2014 - Поле 7 кГ   <br />
2) 2 область:
   1783670 событий в 4 диапазоне заходов <br />
   
   25.09.2014 - Поле 6 кГс  <br />
   старт - RUN=20013    дата: 25.09.2014      ->   cosm_runs_sep14_3.root -  <br />
   финиш - RUN=20017    дата: 25.09.2014 <br />
   старт - RUN=20020    дата: 25.09.2014      ->   cosm_runs_sep14_4.root -   <br />
   финиш - RUN=20023    дата: 25.09.2014 <br />

   1.10.2014 - Поле 6 кГс. <br />

   старт - RUN=20035    дата:  1.10.2014      ->   cosm_runs_oct14_1.root +   <br />
   финиш - RUN=20039    дата:  1.10.2014 <br />

   старт - RUN=20040    дата:  1.10.2014      ->   cosm_runs_oct14_2.root +  <br />
   финиш - RUN=20142    дата:  8.10.2014                        <br />
   8.10.2014 - Поле 6 кГс. <br />
   cosm_runs_oct14_3.root + <br />
   старт - RUN=20206    дата: 12.10.2014      <br />
   финиш - RUN=20211    дата: 12.10.2014      ->   cosm_runs_oct14_4.root + <br />
   cosm_runs_oct14_5.root + <br />
      
      
   1874398 событий в заходах: <br />
       
   старт - RUN=20312    дата: 20.10.2014      ->   cosm_runs_oct14_6.root +  <br />
   cosm_runs_oct14_7.root +  <br />
   старт - RUN=20466    дата: 06.11.2014      ->   cosm_runs_nov14_1.root +  <br /> 
   cosm_runs_nov14_2.root +  <br />
   старт - RUN=20523    дата: 14.11.2014      ->   cosm_runs_nov14_3.root +  <br />
   cosm_runs_nov14_4.root +   <br />  
   старт - RUN=20562    дата: 19.11.2014  <br />
   финиш - RUN=20570    дата: 20.11.2014      ->   cosm_runs_nov14_5.root +    <br />

3) 3 область:  
       1363775 событий в заходах  <br />
      
   финиш - RUN=20602    дата: 26.11.2014      ->   cosm_runs_nov14_6.root +  <br />     
   финиш - RUN=20611    дата: 28.11.2014      ->   cosm_runs_nov14_7.root +  <br />
   cosm_runs_dec14_1.root +  <br />
   cosm_runs_dec14_2.root +   <br />
   RUN=20821    дата: 25.12.2014              ->   cosm_runs_dec14_3.root -- до конца этих заходов поле 6 кГс (лучше смотреть поле по датчику (Холла)  <br />

   1697830 событий в заходах:  <br />

   старт - RUN=20850    дата: 29.12.2014   ->       cosm_runs_dec14_4.root +  <br />
   cosm_runs_dec14_5.root +  <br />
   cosm_runs_jan15_1.root +  <br />
   cosm_runs_jan15_2.root +  <br />
        
   старт - RUN=21101    дата: 02.02.2015   ->       cosm_runs_feb15_1.root +  <br />
   старт - RUN=21115    дата: 03.02.2015   ->       cosm_runs_feb15_2.root +  <br />
                           финиш - RUN=21140    дата: 05.02.2015 <br />
   финиш - RUN=21150    дата: 07.02.2015   ->       cosm_runs_feb15_3.root + <br />
     
   c 08.02.2015 по 12.02.2015 - подъем поля        <br />
   старт - RUN=21194    дата: 13.02.2015   ->       cosm_runs_feb15_4.root +  <br />
      
4) 4 область:        
   1181512 событий в диапазоне заходов  <br />
      
   18.03.2015 по 04.04.2015 - отсутствует поле  <br />
      
   старт - RUN=21548    дата: 05.04.2015   ->       cosm_runs_apr15_1.root  <br />
        
   старт - RUN=21596    дата: 12.04.2015   ->       cosm_runs_apr15_2.root <br />
      
   старт - RUN=21897    дата: 21.05.2015   ->       cosm_runs_may15_1.root    <br />
   ->       cosm_runs_may15_2.root  <br />
      
   c 9.07.2015 - вывод поля  <br />
            
# 2015 год

5) 5 область:         
   2194320 событий в диапазоне заходов  <br />
     
   июнь      1685827 событий в диапазоне заходов  <br />
   cosm_runs_jun15_1.root        -       RUN=22036    -   Поле 6кГс до 09.07.2015    <br />
   cosm_runs_jun15_2.root  <br />
   cosm_runs_jun15_3.root  <br />
   cosm_runs_jun15_4.root  <br />
   cosm_runs_jun15_5.root   <br />          				
   июль       508493 событий в диапазоне заходов     <br /> 
   cosm_runs_jul15_1.root                RUN=22248 -    2015-07-01 -  Поле до 09.07.2015     <br /> 
   cosm_runs_jul15_2.root                RUN=22330  <br />
     
6) 6 область:     
   2052809 событий в диапазоне заходов  <br />

   октябрь                      -         RUN=22370  - 2015-10-01  -     Поле 6кГс с 01.10.2015   <br /> 
   cosm_runs_oct15_1.root  <br />
   cosm_runs_oct15_2.root  <br />
   cosm_runs_oct15_3.root  <br />
   cosm_runs_oct15_4.root  <br />
   cosm_runs_oct15_5.root  <br />
   cosm_runs_nov15_1.root      -          RUN=22576  - 2015-11-01  <br />
   RUN=22628   <br />
 
7) 7 область:   
   1859700 событий в диапазоне заходов  <br />

   cosm_runs_nov15_2.root  <br />
   cosm_runs_nov15_3.root  <br />
   cosm_runs_nov15_4.root  <br />
   cosm_runs_dec15_1.root             -   RUN=22803   -  2015-12-01  <br />
   cosm_runs_dec15_2.root <br />
   cosm_runs_dec15_3.root             -   RUN=22977   =   Поле 6кГс до 24.12.2015  <br />
     
# 2016 год   
   январь                       -         RUN=22996  - 2016-01-15      Поле 6кГс с 15.01.2016   <br /> 
   cosm_runs_jan16_1.root  <br />
   cosm_runs_jan16_2.root   <br />
   cosm_runs_jan16_3.root      -          RUN=23108  <br />

8) 8 область:
   февраль                             -  RUN=23135 - 2016-02-01  <br /> 
   cosm_runs_feb16_1.root  <br />
   cosm_runs_feb16_2.root  <br />
   cosm_runs_feb16_3.root  <br />
   cosm_runs_feb16_4.root  <br />
   cosm_runs_feb16_5.root   <br />    

   1795309 событий в диапазоне заходов    <br />
   cosm_runs_feb16_6.root              RUN=23310  <br />
   cosm_runs_feb16_7.root  <br />
   март                            <br />        
   cosm_runs_march16_1.root  <br />
   cosm_runs_march16_2.root  <br />
   cosm_runs_march16_3.root  <br />
   cosm_runs_march16_5.root       -    RUN=23553    <br />

9) 9 область:
   2134548 событий в диапазоне заходов   <br /> 
   апрель                              -  RUN=23598  <br />
   cosm_runs_apr16_1.root  <br />
   cosm_runs_apr16_2.root  <br />
   cosm_runs_apr16_3.root  <br />
   cosm_runs_apr16_4.root  <br />
   cosm_runs_may16_1.root  <br />
   cosm_runs_may16_2.root        -        RUN=23842  <br />

10) 10 область:      
   2629564 событий в диапазоне заходов   <br />
   май                              -     RUN=23863  <br />

   cosm_runs_may16_3.root  <br />
   cosm_runs_may16_4.root  <br />
   cosm_runs_may16_5.root  <br />
   cosm_runs_may16_6.root  <br />
   июнь                          -        RUN=23998 <br />
   RUN=24079    -  Поле 6кГс до 12.06.2016    <br />
   cosm_runs_jun16_1.root  <br />
   cosm_runs_jun16_2.root  <br />
   cosm_runs_jun16_3.root        -        RUN=24073  <br />

# 2017 год
   RUN=24595   -   Поле 6кГс с 07.02.2017   <br /> 

   RUN=24714  <br />
   cosm_runs_feb17_1.root  -   калориметр все обрезает  <br />
   cosm_runs_feb17_2.root  <br />

   955927 событий  <br />
   cosm_runs_march17_1.root  <br />   
   cosm_runs_march17_2.root      RUN=24925  <br />
   cosm_runs_march17_3.root      RUNstop=25050  <br />
  
   1129600 events  <br />
   cosm_runs_apr17_1.root  <br />
   cosm_runs_apr17_2.root   - RUN=25256       <br />
   cosm_runs_apr17_3.root    - RUNstop=25316  <br />
   
   426035 events  <br />
   cosm_runs_may17_1.root  <br />
   cosm_runs_may17_2.root     -   RUN=25577  <br />
   cosm_runs_may17_3.root     -  RUNstop=25587    <br />
   
   cosm_runs_oct17_1.root -  <br />
   cosm_runs_oct17_2.root -       RUN=25780  <br />
   cosm_runs_oct17_3.root -  <br />
   cosm_runs_oct17_4.root - <br />
   
   cosm_runs_oct17_5.root - <br />
   cosm_runs_oct17_6.root - <br />
   cosm_runs_oct17_7.root - <br />
   cosm_runs_oct17_8.root - <br />
   cosm_runs_oct17_9.root - <br />
   
   cosm_runs_nov17_1.root   -    RUN=26025 <br />
   
11) 11 область:      
   1685846 events <br />
   cosm_runs_dec17_1.root <br /> 
   cosm_runs_dec17_2.root  <br />
   cosm_runs_dec17_3.root     -   RUN=26107 <br />
   cosm_runs_dec17_4.root  <br />
   cosm_runs_dec17_5.root  <br />
   cosm_runs_dec17_6.root  <br />
   cosm_runs_dec17_7.root     -   RUNstop=26320  <br />

   1219009 events <br />
   cosm_runs_jan18_1.root <br /> 
   cosm_runs_jan18_2.root  <br />
   cosm_runs_jan18_3.root      -   RUN=26422 <br />
   cosm_runs_jan18_4.root  <br />
   cosm_runs_jan18_5.root  <br />
   cosm_runs_jan18_6.root  <br />
   cosm_runs_jan18_7.root       -  RUNstop=26611 <br />

12) 12 область:      
   1039271 events <br />
   cosm_runs_feb18_1.root <br /> 
   cosm_runs_feb18_2.root  <br />
   cosm_runs_feb18_3.root    -     RUN=26776 <br />
   cosm_runs_feb18_4.root  <br />
   cosm_runs_feb18_5.root  <br />
   cosm_runs_feb18_6.root  <br />
   cosm_runs_feb18_7.root     -    RUNstop=26949 <br />
 
   921466 events <br />
   cosm_runs_march18_1.root <br /> 
   cosm_runs_march18_2.root  <br />
   cosm_runs_march18_3.root      -   RUN=27085 <br />
   cosm_runs_march18_4.root  <br />
   cosm_runs_march18_5.root      -  RUNstop=27180 <br />

13) 13 область:
   996682 events  <br />   
   cosm_runs_apr18_1.root  <br /> 
   cosm_runs_apr18_2.root  <br />
   cosm_runs_apr18_3.root        -   RUN=27390 <br />
   cosm_runs_apr18_4.root  <br />
   cosm_runs_apr18_5.root        -   RUNstop=27483 <br />

   1618676 events <br />
   cosm_runs_may18_1.root  <br />
   cosm_runs_may18_2.root  <br />
   cosm_runs_may18_3.root        -   RUN=27644 <br />
   cosm_runs_may18_4.root  <br />
   cosm_runs_may18_5.root  <br />
   cosm_runs_may18_6.root        -   RUNstop=27709 <br />

   804907 events <br />
   cosm_runs_jun18_1.root <br /> 
   cosm_runs_jun18_2.root  <br />
   cosm_runs_jun18_3.root        -   RUN=27861 <br />
 
