import pandas as pd
import glob
import numpy as np
import random
from sklearn.model_selection import train_test_split, GridSearchCV, KFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
import sklearn.metrics as metrics



def train_test_split_leaveoneout(x, y, index, stratify=False, random_state=1103):
    idx = list( range(0, len(y)) )
    test_idx = index
    idx.remove(test_idx)
    train_idx = idx

    train_x = x.loc[x.index[train_idx], :]
    train_y = [y[i] for i in train_idx]
    test_x = x.loc[x.index[test_idx], :]
    test_y = y[test_idx]
 
    if stratify:
        random.seed(random_state)
        train_nonresponder_idx = [ i for i, response in enumerate(train_y) if response==0 ]
        train_responder_idx = [ i for i, response in enumerate(train_y) if response==1 ]

        sample_train_nonresponder_idx = sorted( random.sample(train_nonresponder_idx, train_y.count(1)) )
        sample_train_x_nonresponder = train_x.loc[train_x.index[sample_train_nonresponder_idx] , :]
        train_x_responder = train_x.loc[ train_x.index[train_responder_idx], :]
        
        sample_train_y_nonresponder = [ train_y[i] for i in sample_train_nonresponder_idx ]
        train_y_responder = [ train_y[i] for i in train_responder_idx ]
        
        train_x = pd.concat([sample_train_x_nonresponder, train_x_responder], axis=0)
        train_y = sample_train_y_nonresponder + train_y_responder

    return train_x, test_x, train_y, test_y


def Leaveoneout_CV(x, y, feature_name, ML="LogisticRegression"):
    if len(x) != len(y): raise RuntimeError

    pred_prob_list, true_class_list, sample_list, weight_list = [], [], [], []
    for i in range(0, len(y)):

        train_x, test_x, train_y, test_y = train_test_split_leaveoneout(x, y, index=i)

        sample_list.append(test_x.name)
        test_x = np.array(test_x).reshape(1, -1)

        if ML == "LogisticRegression":
            model = LogisticRegression(solver="liblinear", max_iter=100, random_state=1103, class_weight="balanced") # penalty="l1" , class_weight="balanced"

        model.fit(train_x, train_y) 
        

        ## test
        pred_y = model.predict_proba(test_x)[0][1]
        
        pred_prob_list.append(pred_y)
        true_class_list.append( test_y )
        weight_list.append( pd.DataFrame(model.coef_, index=[sample_list[i]], columns=list(train_x.columns)).T )
    AUC = calculate_auc(true_class_list, pred_prob_list)
    AUPRC = metrics.average_precision_score(true_class_list, pred_prob_list)

    feature_weight = pd.concat(weight_list, axis=1)
    pred_perf = pd.DataFrame([AUC, AUPRC], index=["AUC", "AUPRC"], columns=[feature_name]).T
    predicted_result = pd.DataFrame([pred_prob_list, true_class_list, [feature_name] * len(y)], 
                                    index=["pred_prob", "true_class", "feature_name"], columns=sample_list).T
    print("LOOCV completed")

    return pred_perf, predicted_result, feature_weight


def calculate_auc(true_class, pred_prob):
    fpr, tpr, thresholds = metrics.roc_curve(true_class, pred_prob)
    return metrics.auc(fpr, tpr)




