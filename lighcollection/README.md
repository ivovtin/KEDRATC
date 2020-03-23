# Вычисление коэффициентов светосбора для программы моделирования счетчиков АЧС

Для вычисления коэффициентов светосбора данные полученные на космических мюонах разбиваются на временные промежутки 1.5 месяца. <br \>
light_collection - Код для получения коэффициентов светосбора для каждой из области разбиения счетчика (28,58,68 - областей) - усреднение по 5 типам счетчикам. <br \> 
Для запуска на batch: <br \>
```
qsub -shell n -b y -V -N reg1 -cwd ./light_collection 1 1 11022019
```
light_col_coef_cnt - Вычисление поправочного коэффициента для каждого счетчика от среднего значения. <br \>
Для запуска:
``` 
/light_col_coef_cnt 16 1 12122018
```