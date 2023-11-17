# Dependencies
import pandas as pd
import scipy.stats as stats
import numpy as np

# Python built-in dependencies
from itertools import combinations
from itertools import product
from functools import partial
import math


class L2GEQ:
    def __init__(self, csv_path_or_pandas_dataframe, sep=','):
        # 모든 결측값은 0으로 표현한다. 이후 significance test에 np.nan이 나오면 모든 값이 0이거나
        if type(csv_path_or_pandas_dataframe) == str:
            self.data = pd.read_csv(csv_path_or_pandas_dataframe, sep=sep, index_col=0, header=0).fillna(0)
        elif type(csv_path_or_pandas_dataframe) == pd.DataFrame:
            self.data = csv_path_or_pandas_dataframe.fiilna(0)
        else:
            raise ValueError("L2 didn't receive legitimate CSV or DataFrame.")

    # -----------------------------------------------------------------------------
    # House-only functions for L2

    ## For calculation
    def __inverse(p):
        return 1 / p

    def __log10(p):
        return -math.log10(p)

    def __none(p):
        return 1

    def __limit1_linear(p):
        return -p + 1

    def __limit1_circular(p):
        return 1 - math.sqrt(1 - ((p - 1) ** 2))

    def __multiply(p, q):
        return p * q

    def __add(p, q):
        return p + q

    ## Function pointer for statistical significance
    __significance_tests = {
        "ttest": partial(stats.ttest_ind, equal_var=False),
        "mannwhitney": stats.mannwhitneyu,
    }

    ## Function pointer for averaging method
    __averaging_methods = {
        "mean": np.mean,
        "harmonic": stats.hmean,
    }

    ## Function pointer for weighting metric types
    __weighting_methods = {
        "inverse_p": __inverse,
        "log10": __log10,
        "none": __none,
        "limit1_linear": __limit1_linear,
        "limit1_circular": __limit1_circular,
    }

    ## Function pointer for weighting
    __weighting_denominator = {
        "add": __add,
        "multiply": __multiply,
    }

    # For data processing
    ## A built-in function for the naming rule of experiment elements combinations
    def __comb_naming_rule(self, i, j):
        return str(i) + '-' + str(j)

    ## A built-in function for finding indices with the matching label
    def __label_indices(self, i):
        return [index for index in range(len(self.labels_received)) if self.labels_received[index] == i]

    # Split the group in half, returning indices
    def __numbered_indices_combination_without_duplicate(self, the_list):
        cal_list = []
        for i in combinations(the_list, int(len(the_list) / 2)):
            if len(the_list) % int(len(the_list) / 2) == 0:
                cal_list.append(sorted([i, tuple([k for k in the_list if k not in i])]))
            else:
                cal_list.append(sorted([i, tuple([k for k in the_list if k not in i])]))

        else:
            if len(the_list) % int(len(the_list) / 2) == 0:
                cal_list = cal_list[0:int(len(cal_list) / int(len(the_list) / 2))]
        return cal_list

    # Thanks to princekumaras  https://www.geeksforgeeks.org/python-ordered-set/
    def __removeduplicate(self, data):
        countdict = {}
        for element in data:
            if element in countdict.keys():

                # increasing the count if the key(or element)
                # is already in the dictionary
                countdict[element] += 1
            else:
                # inserting the element as key  with count = 1
                countdict[element] = 1
        cal = []
        for key in countdict.keys():
            cal.append(key)
        return cal

    def __sign(self, x):
        if x == 0:
            return 0
        return math.copysign(1, x)

    ## L2
    def __calculate(self, label1, label2, sig_test, avg_method, weight_method, weight_denom="multiply"):
        cal = []
        if sig_test == "mannwhitney":
            cal.extend(
                [
                    self.__weighting_denominator[weight_denom](
                        (self.__averaging_methods[avg_method](self.data.iloc[k, label1]) - self.__averaging_methods[
                            avg_method](self.data.iloc[k, label2])) ** 2,
                        self.__weighting_methods[weight_method](
                            self.__significance_tests[sig_test](self.data.iloc[k, label1],
                                                                self.data.iloc[k, label2]).pvalue)
                    )
                    if (np.isnan(self.data.iloc[k, label1]).any() + np.isnan(self.data.iloc[k, label2]).any() == 0)
                       and (np.sum(self.data.iloc[k, label1] != 0) + np.sum(self.data.iloc[k, label2] != 0) != 0)
                    else np.nan
                    for k in range(self.data.shape[0])
                ]
            )
            return cal
        cal.extend(
            [
                self.__weighting_denominator[weight_denom](
                    (self.__averaging_methods[avg_method](self.data.iloc[k, label1]) - self.__averaging_methods[
                        avg_method](self.data.iloc[k, label2])) ** 2,
                    self.__weighting_methods[weight_method](
                        self.__significance_tests[sig_test](self.data.iloc[k, label1],
                                                            self.data.iloc[k, label2]).pvalue)
                )
                if math.isnan(self.__significance_tests[sig_test](self.data.iloc[k, label1],
                                                                  self.data.iloc[k, label2]).pvalue) == False
                else np.nan
                for k in range(self.data.shape[0])
            ]
        )
        return cal

    def __calculate_with_gsea(self, label1, label2, sig_test, avg_method, weight_method, weight_denom="multiply"):
        cal = []
        cal_gsea = []
        if sig_test == "mannwhitney":
            for k in range(self.data.shape[0]):
                if (np.isnan(self.data.iloc[k, label1]).any() + np.isnan(self.data.iloc[k, label2]).any() == 0) and (
                        np.sum(self.data.iloc[k, label1] != 0) + np.sum(self.data.iloc[k, label2] != 0) != 0):
                    cal.append(
                        self.__weighting_denominator[weight_denom](
                            (self.__averaging_methods[avg_method](self.data.iloc[k, label1]) - self.__averaging_methods[
                                avg_method](self.data.iloc[k, label2])) ** 2,
                            self.__weighting_methods[weight_method](
                                self.__significance_tests[sig_test](self.data.iloc[k, label1],
                                                                    self.data.iloc[k, label2]).pvalue))
                    )
                    cal_gsea.append(self.__sign(
                        self.__averaging_methods[avg_method](self.data.iloc[k, label1]) - self.__averaging_methods[
                            avg_method](self.data.iloc[k, label2])))
                else:
                    cal.append(np.nan)
                    cal_gsea.append(np.nan)
            return cal, cal_gsea

        for k in range(self.data.shape[0]):
            if math.isnan(self.__significance_tests[sig_test](self.data.iloc[k, label1],
                                                              self.data.iloc[k, label2]).pvalue) == False:
                cal.append(
                    self.__weighting_denominator[weight_denom](
                        (self.__averaging_methods[avg_method](self.data.iloc[k, label1]) - self.__averaging_methods[
                            avg_method](self.data.iloc[k, label2])) ** 2,
                        self.__weighting_methods[weight_method](
                            self.__significance_tests[sig_test](self.data.iloc[k, label1],
                                                                self.data.iloc[k, label2]).pvalue))
                )
                cal_gsea.append(self.__sign(
                    self.__averaging_methods[avg_method](self.data.iloc[k, label1]) - self.__averaging_methods[
                        avg_method](self.data.iloc[k, label2])))
            else:
                cal.append(np.nan)
                cal_gsea.append(np.nan)
        return cal, cal_gsea

    # Warning: Assumes you to have Ensembl ID
    def __convert_ensembl_to_hgnc_gene_symbol(self, data, keep_which_duplicate='maximum'):
        # Load the Ensembl id to gene symbol database
        annotation_db = pd.read_csv(self.annotation_db_path, sep=self.__annotation_db_sep)
        annotation_db = annotation_db.rename(columns={"Probe Set ID": "Ensembl_id", "Gene Symbol": "Gene"})[
            ["Ensembl_id", "Gene"]]

        # Convert
        cal = data.merge(annotation_db,
                         how="outer",
                         on="Ensembl_id"
                         )

        if keep_which_duplicate != 'maximum':
            raise NotImplementedError(
                "Multiple Ensembl IDs sometimes map to a single HGNC gene symbol due to alternate splice and etc.\
                The convention is to keep the Ensembl ID with the maximum value, although other methods are possible.\
                In this version of LT2M, duplicate handling is limited to the maximum policy."
            )

        dropped = cal[['Gene', data.keys()[0]]].dropna().sort_values(ascending=False,
                                                                     by=data.keys()[0]).drop_duplicates(
            keep='first')
        final = dropped[['Gene', data.keys()[0]]]
        return final

    # -----------------------------------------------------------------------------

    def designate_label(self, a_list):
        self.labels_received = a_list

    def calculate_L2GEQ(self, sig_test='ttest', avg_method='mean', weight_method='inverse_p', weight_denom='multiply'):
        # Checking group label integrity
        try:
            type(self.labels_received)
        except:
            raise ValueError("L2 didn't receive proper labels. Use designate_label() method.")

        # Broadcast indices to calculate
        cal = {}
        cal_gsea = {}
        for i, j in combinations(self.__removeduplicate(self.labels_received), 2):
            group_a = self.__label_indices(i)
            group_b = self.__label_indices(j)

            cal[self.__comb_naming_rule(i, j)], cal_gsea[self.__comb_naming_rule(i, j)] = self.__calculate_with_gsea(
                group_a, group_b, sig_test, avg_method, weight_method,
                weight_denom)
        else:
            self.weighted_L2 = pd.DataFrame(cal).set_index(self.data.index)
            self.matched_signs_for_L2 = pd.DataFrame(cal_gsea).set_index(self.data.index)

        the_results = {}
        print('---------------', sig_test, avg_method, weight_method, weight_denom, '---------------')
        for i, j in combinations(self.__removeduplicate(self.labels_received), 2):
            cal = self.weighted_L2[self.__comb_naming_rule(i, j)].dropna()
            print(
                i + ' vs. ' + j,
                ":\t",
                math.sqrt(cal.sum() / cal.shape[0])
            )
            the_results[self.__comb_naming_rule(i, j)] = [math.sqrt(cal.sum() / cal.shape[0])]

        self.result = pd.DataFrame(the_results)
        self.parameter_used = sig_test + ' ' + avg_method + ' ' + weight_method + ' ' + weight_denom

    def calculate_withinL2GEQ(self, sig_test='ttest', avg_method='mean', weight_method='inverse_p',
                              weight_denom='multiply'):

        # Checking if the data contains experiments groups that contain +4 elements
        flag = False
        labels_to_iterate = []
        for i in self.__removeduplicate(self.labels_received):
            if len(self.__label_indices(i)) > 3:
                flag = True
                labels_to_iterate.append(i)
        if flag == False:
            return "No groups contained more than 3 members, making it impossible for each group of combination to contain 2+ members."

        within_L2GEQ = {}

        for a in labels_to_iterate:
            for b in self.__numbered_indices_combination_without_duplicate(self.__label_indices(a)):
                cal = []
                cal.extend(self.__calculate(list(b[0]), list(b[1]), sig_test, avg_method, weight_method, weight_denom))
                within_L2GEQ[self.__comb_naming_rule(a, b)] = cal.copy()
            else:
                self.inner_L2GEQ = pd.DataFrame(within_L2GEQ)

        print('---------------', sig_test, avg_method, weight_method, weight_denom, '---------------')
        for i in self.inner_L2GEQ:
            cal = self.inner_L2GEQ[i].dropna()
            print(i,
                  '\t',
                  math.sqrt(cal.sum() / cal.shape[0])
                  )

    def calculate_L2GEQ_custom_index(self, custom_list1, custom_list2, sig_test='ttest', avg_method='mean',
                                     weight_method='inverse_p', weight_denom='multiply'):
        # Checking group label integrity
        try:
            type(self.labels_received)
        except:
            raise ValueError("L2 didn't receive proper labels. Use designate_label() method first.")

        cal = pd.DataFrame(
            self.__calculate(custom_list1, custom_list2, sig_test, avg_method, weight_method, weight_denom)).dropna()

        print(math.sqrt(cal.sum() / cal.shape[0]))

    def parameter_iterator(self, mode='main', sig_test=['ttest', 'mannwhitney'], avg_method=['mean', 'harmonic'],
                           weight_method=['inverse_p', 'log10', 'none', 'limit1_linear', 'limit1_circular'],
                           weight_denom=['multiply'], extract_degs_too=True):
        print(mode)
        mode_dict = {'main': self.calculate_L2GEQ, 'within': self.calculate_withinL2GEQ}
        for i in product(sig_test, avg_method, weight_method, weight_denom, repeat=1):
            mode_dict[mode](i[0], i[1], i[2], i[3])
            if mode == 'main' and extract_degs_too == True:
                self.extract_DEGs()

    def set_annotation_db(self, path, sep):
        self.annotation_db_path = path
        self.__annotation_db_isGiven = True
        self.__annotation_db_sep = sep

    def extract_DEGs(self, keys='all', cutoff=0.0005, cutoff_type='relative', save_file=True, verbose=True,
                     path_annotation_db=None, annotation_db_sep='\t'):
        try:
            type(self.weighted_L2)
        except:
            raise AttributeError("You should calculate L2 first to extract DEGs.")

        annotation_flag = False

        if path_annotation_db == None and self.__annotation_db_isGiven == True:
            annotation_flag = True
        elif path_annotation_db != None:
            annotation_flag = True
            self.annotation_db_path = path_annotation_db
            self.__annotation_db_isGiven = True
            self.__annotation_db_sep = annotation_db_sep

        if keys == 'all':
            keys = self.weighted_L2.keys()
        else:
            try:
                if type(keys) != list:
                    raise ValueError("Argument keys should be a list.")
                if np.isin(keys, self.weighted_L2.keys()).sum() == len(keys):
                    cal = self.weighted_L2[keys]
                else:
                    raise ValueError("One or more of your keys are not found in the data.")
            except:
                raise AttributeError("There is an error in your input to argument keys.")

        cal_running_sum = {}
        for i in keys:
            cal_length = self.weighted_L2[i].dropna().shape[0]
            cal_running_sum[i] = self.weighted_L2[i].dropna().sort_values(ascending=False).cumsum().apply(
                lambda x: x / cal_length).apply(math.sqrt)

        self.running_sum = cal_running_sum

        self.DEGs = {}
        if annotation_flag == True:
            self.DEGs_gene_symbol = {}
        for i in keys:
            cal = pd.DataFrame(self.running_sum[i]).dropna().diff()
            cal.iloc[0, 0] = self.running_sum[i][0]
            # GSEA doesn't use cutoff
            gsea_cal = (cal * pd.DataFrame(self.matched_signs_for_L2[i][cal.index])).sort_values(ascending=False, by=i)
            if cutoff_type == 'relative':
                # cal = cal[cal > self.result.iloc[0, :].mean() * cutoff].dropna()
                cal = cal[cal > self.result[i][0] * cutoff].dropna()
            if cutoff_type == 'absolute':
                cal = cal[cal > cutoff].dropna()
            self.DEGs[i] = cal.index

            if annotation_flag == True:
                converted = self.__convert_ensembl_to_hgnc_gene_symbol(pd.DataFrame(cal))
                self.DEGs_gene_symbol[i] = converted["Gene"]

            if save_file == True:
                if annotation_flag == True:
                    converted.to_csv(i + ' ' + self.parameter_used + '.tsv', sep='\t', index=False)
                    self.__convert_ensembl_to_hgnc_gene_symbol(gsea_cal).to_csv("GSEA " + i + ' ' + self.parameter_used + '.tsv', sep='\t', index=False)
                else:
                    pd.DataFrame(cal).to_csv(i + ' ' + self.parameter_used + '.tsv', sep='\t', index=False)
                    gsea_cal.to_csv("GSEA " + i + ' ' + self.parameter_used + '.tsv', sep='\t', index=False)

        if verbose == True:
            if annotation_flag == True:
                for i in self.DEGs_gene_symbol:
                    print("----------", "Differentially expressed genes for :", i, ",", len(self.DEGs_gene_symbol[i]),
                          "genes identified",
                          "----------")
                    for j in self.DEGs_gene_symbol[i]:
                        print(j)
            else:
                for i in self.DEGs:
                    print("----------", "Differentially expressed genes for :", i, ",", len(self.DEGs[i]),
                          "genes identified",
                          "----------")
                    for j in self.DEGs[i]:
                        print(j)