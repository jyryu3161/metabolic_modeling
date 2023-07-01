import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon, ttest_ind, ttest_rel, ranksums#
from scipy.stats import shapiro
from statsmodels.stats.multitest import fdrcorrection


def two_grouped_data_comparison(condition1_df, condition2_df, related=False, p_value_cutoff=0.05, filter_by='FDR',show_missing_data=False):
    # Filter_by
    while filter_by not in ['FDR', 'Bonferroni', 'None']:
        print('Input argument filter_by is not correct.')
        filter_by=input('Choose[\'FDR\'/\'Bonferroni\'/\'None\']:')

    # Filtration of missing data and statistically unsignificant data
        #remove nan
    pure_condition1_df = condition1_df.dropna(how='any',axis=0)
    pure_condition2_df = condition2_df.dropna(how='any',axis=0)
    pure_index_list = list(set(pure_condition1_df.index) & set(pure_condition2_df.index))
        #remove all-equal sample set
    drop_index = []
    for ind in pure_index_list:
        sample_set = list(pure_condition1_df.loc[ind].values)+list(pure_condition2_df.loc[ind].values)
        if len(sample_set) == sample_set.count(min(sample_set)):
            drop_index.append(ind)
    pure_condition1_df = pure_condition1_df.drop(index=drop_index)
    pure_condition2_df = pure_condition2_df.drop(index=drop_index)
    pure_index_list = list(set(pure_condition1_df.index) & set(pure_condition2_df.index))

    missed_df1 = pd.DataFrame(columns = condition1_df.columns,
                                index=list(set(condition1_df.index)-set(pure_index_list)),
                                data=condition1_df)
    missed_df2 = pd.DataFrame(columns = condition2_df.columns,
                                index=list(set(condition2_df.index)-set(pure_index_list)),
                                data=condition2_df)
    # identifying comparison test
    if len(pure_condition1_df.columns) < 10 or len(pure_condition2_df.columns) < 10:
        test_method = 'non_parametric'
    else:
        test_method = 'parametric'
        for each_index in pure_index_list:
            pop1 = list(pure_condition1_df.loc[each_index,:].values)
            pop2 = list(pure_condition2_df.loc[each_index,:].values)
            statistics1, shapiro_p_value1 = shapiro(pop1)
            statistics2, shapiro_p_value2 = shapiro(pop2)
            if shapiro_p_value1 < 0.05 or shapiro_p_value2 < 0.05:
                test_method='non_parametric'
                break

    if test_method == 'non_parametric':
        if related:
            method='Wilcoxon signed-rank test'
            target_function = wilcoxon
        else:
            method='Wilcoxon rank-sum test'
            target_function = ranksums
    elif test_method == 'parametric':
        if related:
            method = 'T-test on two paired samples'
            target_function = ttest_rel
        else:
            method = 'T-test on two independent samples'
            target_function = ttest_ind
    # Statistical comparison
    result_df_dict = {}
    result_df_columns = ['Condition1 mean','Condition1 std','Condition2 mean','Condition2 std','Method','P','log2 (condition2/condition1)','Status']
    result_df_index = []
    for each_col in result_df_columns:
        result_df_dict[each_col] = []

    for each_index in pure_index_list:
        error=False
        pop1 = list(pure_condition1_df.loc[each_index,:].values)
        pop2 = list(pure_condition2_df.loc[each_index,:].values)
        try:
            statistics, pvalue = target_function(pop1, pop2)
        except:
            error=True
        ## in case p value has error
        if pvalue!=pvalue or error==True:
            missed_df1 = pd.concat([missed_df1, pd.DataFrame(pop1, columns=[each_index],index=list(pure_condition1_df.columns.values)).T],axis=0)
            missed_df2 = pd.concat([missed_df2, pd.DataFrame(pop2, columns=[each_index],index=list(pure_condition2_df.columns.values)).T],axis=0)
            pure_condition1_df.drop(index=[each_index],inplace=True)
            pure_condition2_df.drop(index=[each_index],inplace=True)
            continue
        ## Status
        mean1 = np.mean(pop1)
        mean2 = np.mean(pop2)
        if mean1 * mean2 <0.0:
            Status = 'Reversed'
        elif mean1 == mean2:
            Status = 'Unchanged'
        elif abs(mean2) > abs(mean1):
            Status = 'UP'
        else:
            Status = 'DOWN'

        statistical_property = [mean1, np.std(pop1), mean2, np.std(pop2), method, pvalue, np.log2(mean2/mean1), Status]
        result_df_index.append(each_index)
        for each_col in result_df_columns:
            result_df_dict[each_col].append(statistical_property[result_df_columns.index(each_col)])

    # p-value correction
        # FDR correction
    rejected, p_value_corrected = fdrcorrection(result_df_dict['P'])
    result_df_columns.insert(6, 'FDR corrected P')
    result_df_dict['FDR corrected P'] = list(p_value_corrected)
        # Bonferroni
    if result_df_index == 0:
        alpha = p_value_cutoff
    alpha = p_value_cutoff/len(result_df_index)
    result_df_columns.insert(7,'Bonferroni (alpha=%s)'%str(alpha))
    bonferroni_passed = np.array(result_df_dict['P']) < alpha
    bonferroni_result = [str(i).replace('True','Passed').replace('False','Failed') for i in bonferroni_passed]
    result_df_dict['Bonferroni (alpha=%s)'%str(alpha)] = bonferroni_result
    # result_filtration
    if filter_by != 'None':
        show_missing_data=False
        if filter_by=='FDR':
            filter_array = np.array(result_df_dict['FDR corrected P'])<p_value_cutoff
        elif filter_by=='Bonferroni':
            filter_array = np.array(result_df_dict['Bonferroni (alpha=%s)'%str(alpha)]) == 'Passed'
        for each_col in result_df_columns:
            result_df_dict[each_col] = list(np.array(result_df_dict[each_col])[filter_array])
        result_df_index = list(np.array(result_df_index)[filter_array])

    if not show_missing_data:
        result_df = pd.DataFrame(data = result_df_dict, columns = result_df_columns, index=result_df_index)
        return result_df

    # Addition of missing data
    print('Adding missing data...')
    missed_index = list(set(missed_df1.index.values).union(set(missed_df2.index.values)))
    for each_index in missed_index:
        statistical_property = [float('nan')]*len(result_df_dict.keys())
        pop1 = []
        pop2 = []
        ## filter nan data
        if each_index in missed_df1.index:
            raw_pop1 = list(missed_df1.loc[each_index].values)
            for each_value in raw_pop1:
                if each_value == each_value:
                    pop1.append(each_value)
        if each_index in missed_df2.index:
            raw_pop2 = list(missed_df2.loc[each_index].values)
            for each_value in raw_pop2:
                if each_value == each_value:
                    pop2.append(each_value)
        if pop1:
            statistical_property[0:2] = [np.mean(pop1), np.std(pop1)]
        if pop2:
            statistical_property[2:4] = [np.mean(pop2), np.std(pop2)]
        if pop1 and pop2:
            statistical_property[8] = np.log2(np.mean(pop2)/np.mean(pop1))
        for each_col in result_df_columns:
            result_df_dict[each_col].append(statistical_property[result_df_columns.index(each_col)])
        result_df_index.append(each_index)

    result_df = pd.DataFrame(data = result_df_dict, columns = result_df_columns, index=result_df_index)
    return result_df

