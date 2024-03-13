
# Identifying the first Open Reading Frame

In this program, we are identifying the first orf and translating first identified orf of given fasta file 


**Created by: Aruna Kumar**

**Created Date: Feb-23-2023**


```
``$python3  translate_first_orf.py data/dmel-all-chromosome-r6.17.fasta

```


## 1. Finding a first orf of [dmel-all-chromosome-r6.17.fasta]
```shell
(venv) aruna@Pradeeps-MacBook-Pro module-04-aruna-kumar % python3 translate_first_orf.py data/dmel-all-chromosome-r6.17.fasta

```
# output
```
2L:      MHDRGSRTDI*
2R:      MSF*
3L:      MIAYARVVPTYCAL*
3R:      MGPSTEPSTEPSTGPVRDQYGTSTGPVRDQYGTSTGPVRDQYGTSTGPVRDQYGTSTGPSTEPSTGPVRDQYGTSTGPVRD*
4:       MNGIIIGNSTIFNNLYQSTNLVNALLIYLT*
```


## Test behavior of translate_first_orf.py ##

```
(venv) aruna@Pradeeps-MacBook-Pro module-04-aruna-kumar % pytest test_translate_first_orf.py   
```
# output
 
```
 ========(venv) aruna@Pradeeps-MacBook-Pro module-04-aruna-kumar % pytest test_translate_first_orf.py
======================================================================== test session starts ========================================================================
platform darwin -- Python 3.8.9, pytest-7.2.1, pluggy-1.0.0
rootdir: /Users/aruna/PyCharm/module-04-aruna-kumar
collected 5 items                                                                                                                                                   

test_translate_first_orf.py .....                                                                                                                             [100%]

========================================================================= 5 passed in 0.41s ================================
```




