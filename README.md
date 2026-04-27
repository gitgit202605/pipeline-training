# pipeline-training
 how to solve prevalence question

#1. Simple count point estimate, using 600 IHC data as ground truth
### calculate the raw percentage of both cutoffs.

# cutoff1: 2+ 50%, count samples where intensity >= 2+ AND coverage >= 50%
# cutoff2: 1+ 25%, count samples where intensity >= 1+ AND coverage >= 25%

##Observed_Prevalence = Number of Samples meeting Cutoff / Total samples *100


#2. The statistical inference confidence interval

from statsmodels.stats.proportion import proportion_confint
from scipy.stats import binomtest


total_n = len(joint_df)
pos_cutoff_a = len(joint_df[(joint_df['tumor_cell_intensity'] >= 2) & (joint_df['tumor_cell_percent'] >= 50)])
pos_cutoff_b = len(joint_df[(joint_df['tumor_cell_intensity'] >= 1) & (joint_df['tumor_cell_percent'] >= 25)])


results = []

for label, count, definition in [('Cutoff A',pos_cutoff_a, '2+ in 50%' ), ('Cutoff B',pos_cutoff_b, '1+ in 25%')]:

    # use scipy.stats clopper pearson beta method
    #prevalence = count / total_n 
    #ci_low, ci_upp = proportion_confint(count,total_n,alpha = 0.05,method='beta')
    #print(f"{label} prevelance : {prevalence:.1%} (95% CI: [{ci_low:.1%} - {ci_upp:.1%}])")

    # use scipy.stats.binomtest for the exact method
    bt = binomtest ( k = count, n = total_n)
    ci = bt.proportion_ci(confidence_level = 0.95, method = 'exact')

    results.append({
        'Cutoff' : label,
        'Definition' : definition,
        'Positives' : count,
        'Prevalence' : bt.statistic,
        'CI_Low' : ci.low,
        'CI_High' : ci.high      
    })
    
results_df = pd.DataFrame(results)

#3 visualize, side by side 
sns.set_theme(style='whitegrid')
fig,(ax1,ax2) = plt.subplots(1,2,figsize = (16,6), gridspec_kw={'width_ratios': [1,1.5]})

#plot A: point plot with Error Bars
ax1.errorbar(results_df['Cutoff'],results_df['Prevalence'],
             yerr = [results_df['Prevalence'] - results_df['CI_Low'],
                     results_df['CI_High'] - results_df['Prevalence']],
             fmt='o',capsize=8,markersize=12,color='#2c7fb8',elinewidth=3)

ax1.set_title(f'CLDN18.2 Prevalence Comparison (N={total_n})',fontsize=18, fontweight='bold')
ax1.set_ylabel('Prevalence Rate',fontsize=18)
ax1.set_ylim(0,results_df['CI_High'].max() + .1)

# add text labels for percentages

for i, row in results_df.iterrows():
    ax1.text(i, row['Prevalence'] + 0.03, f"{row['Prevalence']:.1%}", ha='center', fontweight='bold')

#plot B: summary table

ax2.axis('off')
table_data = results_df[['Cutoff','Definition','Positives','Prevalence','CI_Low','CI_High']].copy()
#Format for display
table_data['Prevalence'] = table_data['Prevalence'].map('{:.1%}'.format)
table_data['95% CI(Exact)'] = table_data.apply(lambda x: f"[{x['CI_Low']:.1%},{x['CI_High']:.1%}]",axis=1)
table_data['Prevalence'] = table_data.drop(columns = ['CI_Low','CI_High'])

the_table = ax2.table(cellText=table_data.values, colLabels = table_data.columns,
                      loc='center',cellLoc='center')
the_table.auto_set_font_size(False)
the_table.set_fontsize(11)
the_table.scale(1.2,2.5)

plt.tight_layout()
plt.show()

print(table_data.to_string(index=False))


---------------------------------------------------------------------------
ValueError                                Traceback (most recent call last)
/tmp/ipykernel_2081687/214817855.py in ?()
     69 table_data = results_df[['Cutoff','Definition','Positives','Prevalence','CI_Low','CI_High']].copy()
     70 #Format for display
     71 table_data['Prevalence'] = table_data['Prevalence'].map('{:.1%}'.format)
     72 table_data['95% CI(Exact)'] = table_data.apply(lambda x: f"[{x['CI_Low']:.1%},{x['CI_High']:.1%}]",axis=1)
---> 73 table_data['Prevalence'] = table_data.drop(columns = ['CI_Low','CI_High'])
     74 
     75 the_table = ax2.table(cellText=table_data.values, colLabels = table_data.columns,
     76                       loc='center',cellLoc='center')

~/miniconda3/envs/spacec/lib/python3.10/site-packages/pandas/core/frame.py in ?(self, key, value)
   4308             self._setitem_frame(key, value)
   4309         elif isinstance(key, (Series, np.ndarray, list, Index)):
   4310             self._setitem_array(key, value)
   4311         elif isinstance(value, DataFrame):
-> 4312             self._set_item_frame_value(key, value)
   4313         elif (
   4314             is_list_like(value)
   4315             and not self.columns.is_unique

~/miniconda3/envs/spacec/lib/python3.10/site-packages/pandas/core/frame.py in ?(self, key, value)
   4436             loc = self.columns.get_loc(key)
   4437             cols = self.columns[loc]
   4438             len_cols = 1 if is_scalar(cols) or isinstance(cols, tuple) else len(cols)
   4439             if len_cols != len(value.columns):
-> 4440                 raise ValueError("Columns must be same length as key")
   4441 
   4442             # align right-hand-side columns if self.columns
   4443             # is multi-index and self[key] is a sub-frame

ValueError: Columns must be same length as key




