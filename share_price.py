# -*- coding: utf-8 -*-
"""
"""

# -*- coding: utf-8 -*-
"""
LSM 1 - gareth 
"""

import matplotlib . pyplot as plt
import numpy as np
from matplotlib import cm
import math
from scipy import ndimage, misc
from scipy import signal

no_shares = 10000;
buy_price = 2.51;
exchange_rate = 1.1;
max_sell=4;
min_sell=2;
tax=0.52;
cgt=0.33;


N=100;

profit_sell = np.zeros((N,1));
profit_hold = np.zeros((N,1));
price = np.zeros((N,1));

for i in range(0,N):
    # get share price over the interval
    current_val = min_sell + (i)*(max_sell-min_sell)/N;
    price[i][0]=current_val;

    # straight forward profit from selling shares and paying tax
    profit_sell[i][0] = no_shares*(1-tax)*( current_val*exchange_rate );
    
    # get profit after CGT from selling 1 share after purchase at buy_price
    pro = (cgt)*(current_val-buy_price)*exchange_rate ;
    if(pro<0):
        pro=0;
        
    # work out profit from selling shares and adding in the CGT profit
    profit_hold[i][0] = no_shares*(1-tax)*( current_val*exchange_rate - pro);
   
    
fig = plt.figure(figsize=(8, 6)) ;
plt.plot(price,profit_sell);
plt.plot(price,profit_hold);
plt.show();
 