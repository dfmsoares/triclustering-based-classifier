
if __name__ == "__main__":
    import os
    import identify_triclusters as idt
    import sys
    import os
    import numpy as np
    import pandas as pd
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import cross_validate, RepeatedStratifiedKFold, cross_val_score, train_test_split
    from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, make_scorer

    path = sys.argv[1]
    directory = os.fsencode(path)
    target_var = sys.argv[2]
    n = int(sys.argv[3])
    outfile = sys.argv[4]
    n_splits = int(sys.argv[5])
    n_repts = int(sys.argv[6])

    

    table = {"File":[], "AUC": [], "AUC_Std":[], "CA":[], "CA_Std":[], "Sens": [], "Sens_Std":[], "Spec": [], "Spec_Std":[]}
    for file in sorted(os.listdir(directory)):
        filename = os.fsdecode(file)
        if filename.endswith(".csv"): 

            cts = pd.read_csv(path + '/' + filename)

            X = cts.iloc[:,1:len(cts.columns)-1].values
            y = cts[target_var].values

            model = RandomForestClassifier(n_estimators=100, max_depth=7500, random_state=0)

            def tn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 0]
            def fp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 1]
            def fn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[1, 0]
            def tp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[1, 1]

            def sens(y_true, y_pred): return tp(y_true, y_pred)/(fn(y_true, y_pred) + tp(y_true, y_pred))
            def spec(y_true, y_pred): return tn(y_true, y_pred)/(fp(y_true, y_pred) + tn(y_true, y_pred))

            scoring = {'accuracy': make_scorer(accuracy_score), 'roc_auc': 'roc_auc', 'sens': make_scorer(sens), 'spec': make_scorer(spec)}

            rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repts,random_state=36851234)
            scores = cross_validate(model, X, y, scoring=scoring, cv=rskf)

            print ("#"*10)
            print()
            print(" =======", str(n) + "TPS", "=======")
            print("File: " + filename)
            table["File"].append(filename)

            print("AUC: ", scores['test_roc_auc'].mean(), "+-", scores['test_roc_auc'].std())
            table["AUC"].append(scores['test_roc_auc'].mean())
            table["AUC_Std"].append(scores['test_roc_auc'].std())

            print("CA: ", scores['test_accuracy'].mean(), "+-", scores['test_accuracy'].std())
            table["CA"].append(scores['test_accuracy'].mean())
            table["CA_Std"].append(scores['test_accuracy'].std())

            print("Sensitivity: ", scores['test_sens'].mean(), "+-", scores['test_sens'].std())
            table["Sens"].append(scores['test_sens'].mean())
            table["Sens_Std"].append(scores['test_sens'].std())

            print("Specificity: ", scores['test_spec'].mean(), "+-", scores['test_spec'].std())
            table["Spec"].append(scores['test_spec'].mean())
            table["Spec_Std"].append(scores['test_spec'].std())

            print("Feature Importance:")
            model.fit(X,y)
            importances = model.feature_importances_
            indices = np.argsort(importances)[::-1]
            # Print the feature ranking (top20)
            final = 20
            if len(indices) < 20:
                final = len(indices)
            for f in range(final):
                print("%d. feature %s (%f)" % (f + 1, cts.columns[1:len(cts.columns)-1][indices[f]], importances[indices[f]]))

    
    pd.DataFrame(table).to_csv(outfile)
                    