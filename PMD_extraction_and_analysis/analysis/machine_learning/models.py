from config import *

def select_model(y):
    if MODEL == 'rf': 
        from sklearn.ensemble import RandomForestClassifier
        model = RandomForestClassifier(criterion='gini', class_weight='balanced_subsample',  random_state=RANDOM_SEED)

    if MODEL == 'svc':
        from sklearn.svm import SVC
        model = SVC(kernel='linear', probability=True, class_weight='balanced', random_state=RANDOM_SEED)

    if MODEL == 'lr':
        from sklearn.linear_model import LogisticRegression
        model = LogisticRegression(solver='lbfgs', penalty='l2', class_weight='balanced', random_state=RANDOM_SEED, max_iter=5000) 
       
    if MODEL == 'xgb-tree': 
        import xgboost as xgb
        # If we are using LOO validation, we can approximate the scale_pos_weight by the ratio of the number of negative and positive examples
        # Assuming you have the number of positive and negative examples
        num_negative = sum(y == 0) 
        num_positive = sum(y == 1) 
        scale_pos_weight = num_negative / num_positive
        model = xgb.XGBClassifier(booster='gbtree', objective="binary:logistic", scale_pos_weight=scale_pos_weight, random_state=RANDOM_SEED)

    if MODEL == 'xgb-linear':
        import xgboost as xgb 
        num_negative = sum(y == 0)  
        num_positive = sum(y == 1)  
        scale_pos_weight = num_negative / num_positive
        model = xgb.XGBClassifier(booster='gblinear', objective="binary:logistic", scale_pos_weight=scale_pos_weight, random_state=RANDOM_SEED)

    if MODEL == 'knn':
        from sklearn.neighbors import KNeighborsClassifier
        model = KNeighborsClassifier()	

    if MODEL == 'tabpfn':
        from tabpfn import TabPFNClassifier
        model = TabPFNClassifier(n_estimators=25, balance_probabilities=True, random_state=RANDOM_SEED)

    if MODEL == 'phe_tabpfn':
        from tabpfn_extensions.post_hoc_ensembles.sklearn_interface import AutoTabPFNClassifier
        # phe means post hoc ensemble. (post hoc is latin for after the fact).
        # The AutoTabPFNClassifier automatically tries to boost performance by ensembling 100 TabPFN models post-training.
        # You can control the runtime using ´max_time´ and need to make no further adjustments to get best results.
        # but training 100 models for 150 folds takes a lot of time.
        # The results of this are not better than the results of the normal tabpfn model. (at least for max_time=30)
        # I might try this one more time overnight with maxtime = 200,
        # HOW DO WE BALANCE THE CLASSES HERE?
        model = AutoTabPFNClassifier(random_state=RANDOM_SEED, max_time=30)

    if MODEL == 'hpo_tabpfn':
        from tabpfn_extensions.hpo import TunedTabPFNClassifier
        # hpo means hyperparameter optimization. 
        # HOW DO WE BALANCE THE CLASSES HERE?
        model = TunedTabPFNClassifier(n_trials=50, metric='roc_auc', random_state=RANDOM_SEED, verbose=False)

    if MODEL == 'lgbmc':
        import lightgbm as lgb
        num_negative = sum(y == 0)  
        num_positive = sum(y == 1)  
        scale_pos_weight = num_negative / num_positive
        model = lgb.LGBMClassifier(objective='binary', scale_pos_weight=scale_pos_weight, random_state=RANDOM_SEED)

    if MODEL == 'ensemble':
        from sklearn.ensemble import VotingClassifier
        num_negative = sum(y == 0)  
        num_positive = sum(y == 1)  
        scale_pos_weight = num_negative / num_positive
    
        LR_model = LogisticRegression(solver='lbfgs', penalty='l2', class_weight='balanced', random_state=RANDOM_SEED)
        SVC_model = SVC(kernel='linear', probability=True, class_weight='balanced', random_state=RANDOM_SEED)
        XGB_TREE_model = xgb.XGBClassifier(booster='gbtree', objective="binary:logistic", scale_pos_weight=scale_pos_weight, random_state=RANDOM_SEED)
        XGB_LINEAR_model = xgb.XGBClassifier(booster='gblinear', objective="binary:logistic", scale_pos_weight=scale_pos_weight, random_state=RANDOM_SEED)
        model = VotingClassifier(estimators=[('lr', LR_model), ('svc', SVC_model), ('xgb_tree', XGB_TREE_model), ('xgb_linear', XGB_LINEAR_model)], voting='soft')

    return model