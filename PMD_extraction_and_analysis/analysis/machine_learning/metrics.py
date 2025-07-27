import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, roc_auc_score, confusion_matrix, roc_curve, average_precision_score
from config import *
from matplotlib import gridspec
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, confusion_matrix
from sklearn.preprocessing import label_binarize


def metrics_for_report(test_true, test_prob):
    """
    Plot AUROC, AUPRC, confusion matrix and a table with classification metrics side by side in subplots.  
    These were designed for a binary classifier.
    """
    
    # Take probabilities for the positive class
    test_prob = test_prob[:, 1]  

    # Calculate ROC metrics
    roc_auc = roc_auc_score(test_true, test_prob)
    fpr, tpr, _ = roc_curve(test_true, test_prob)

    # Calculate PRC metrics
    precision, recall, _ = precision_recall_curve(test_true, test_prob)
    prc_auc = average_precision_score(test_true, test_prob)

    # Calculate confusion matrix metrics
    cm = confusion_matrix(test_true, (test_prob > 0.5).astype(int))
    tp = cm[1, 1]
    fp = cm[0, 1]
    fn = cm[1, 0]
    tn = cm[0, 0]
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    precision_metric = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    f1 = 2 * (precision_metric * sensitivity) / (precision_metric + sensitivity) if (precision_metric + sensitivity) > 0 else 0

    # Create subplots with gridspec
    # Don't change these sizes, else the fontsize will be wrong
    fig = plt.figure(figsize=(33/2.54, 10/2.54))  # Figure size in inches, I pass the size in cm and convert to inches
    gs = gridspec.GridSpec(1, 3) 
  
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    gs_right = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[2], height_ratios=[1, 1])  # Ensure equal height for c and d
    ax3_top = fig.add_subplot(gs_right[0])
    ax3_bottom = fig.add_subplot(gs_right[1])

    # Plot ROC curve
    ax1.plot(fpr, tpr, label=f'AUC = {roc_auc:.3f}')
    ax1.set_xlabel('False Positive Rate', fontsize=14)
    ax1.set_ylabel('True Positive Rate', fontsize=14)
    ax1.set_title('ROC Curve', fontsize=14)
    ax1.legend(loc='lower right')
    ax1.text(-0.1, 1.15, 'a', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

    # Plot PRC curve
    ax2.plot(recall, precision, label=f'AUC = {prc_auc:.3f}')
    ax2.set_xlabel('Recall', fontsize=14)
    ax2.set_ylabel('Precision', fontsize=14)
    ax2.set_title('PRC Curve', fontsize=14)
    ax2.legend(loc='lower left')
    ax2.text(-0.1, 1.15, 'b', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')
    
    # Top table: Confusion matrix
    cm_cell_text = [[f"{cm[i, j]}" for j in range(cm.shape[1])] for i in range(cm.shape[0])]
    cm_cell_text.insert(0, [''] + ['Pred ' + label for label in ['Stroke', 'Seizure']])  # Add column headers
    for i in range(1, len(cm_cell_text)):
        cm_cell_text[i].insert(0, 'True ' + ['Stroke', 'Seizure'][i - 1])  # Add row headers

    ax3_top.axis('off')
    cm_table = ax3_top.table(
        cellText=cm_cell_text,
        loc='center',
        cellLoc='center'
    )
    cm_table.auto_set_font_size(False)
    cm_table.set_fontsize(14)
    cm_table.scale(1.4, 1.6)  # Adjusted height scaling for better text fitting

    # Make row and column labels bold
    for (row, col), cell in cm_table.get_celld().items():
        if row == 0 or col == 0:  # Header cells
            cell.set_text_props(weight='bold')

    ax3_top.text(-0.2, 1.05, 'c', transform=ax3_top.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

    # Bottom table: Classification metrics
    metrics = ['Specificity', 'Sensitivity', 'Precision', 'NPV', 'F1']
    values = [specificity, sensitivity, precision_metric, npv, f1]
    metrics_cell_text = [[metrics[i], f"{values[i]:.3f}"] for i in range(len(metrics))]

    ax3_bottom.axis('off')
    metrics_table = ax3_bottom.table(cellText=metrics_cell_text, colLabels=["Metric", "Score"], loc='center', cellLoc='center')
    metrics_table.auto_set_font_size(False)
    metrics_table.set_fontsize(14)
    metrics_table.scale(1.2, 1.6)  # Adjusted scaling for better text fitting

    # Make column headers bold
    for (row, col), cell in metrics_table.get_celld().items():
        if row == 0:  # Header row
            cell.set_text_props(weight='bold')

    ax3_bottom.text(-0.2, 1.05, 'd', transform=ax3_bottom.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.3)
    plt.savefig(f'analysis/ml_results.png', dpi=600, bbox_inches='tight')
    plt.show()

    print(" ")
    print(f"AUROC = {roc_auc}")
    print(f"AUPRC = {prc_auc}")
    print(f"Sensitivity = {sensitivity:.4f}")
    print(f"Specificity = {specificity:.4f}")
    print(f"Precision = {precision_metric:.4f}")
    print(f"NPV = {npv:.4f}")
    print(f"F1 = {f1:.4f}")


def metrics_with_ci(test_true, test_prob, nr_bootstraps=5000, interval=0.95):
    """
    Plot AUROC, AUPRC, confusion matrix, and a table with classification metrics side by side in subplots,
    including confidence intervals for AUROC, AUPRC, and other metrics using bootstrapping.
    """
    
    # Take probabilities for the positive class
    test_prob = test_prob[:, 1]  

    # Initialize lists to store bootstrap results
    roc_auc_scores = []
    prc_auc_scores = []
    sensitivity_scores = []
    specificity_scores = []
    precision_scores = []
    npv_scores = []
    f1_scores = []

    # Bootstrap sampling
    for _ in range(nr_bootstraps):
        # Resample the data with replacement
        indices = np.random.choice(len(test_true), size=len(test_true), replace=True)
        y_true_resampled = test_true[indices]
        y_prob_resampled = test_prob[indices]

        # Calculate metrics for the resampled data
        roc_auc = roc_auc_score(y_true_resampled, y_prob_resampled)
        prc_auc = average_precision_score(y_true_resampled, y_prob_resampled)
        cm = confusion_matrix(y_true_resampled, (y_prob_resampled > 0.5).astype(int))
        tp = cm[1, 1]
        fp = cm[0, 1]
        fn = cm[1, 0]
        tn = cm[0, 0]
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        precision_metric = tp / (tp + fp) if (tp + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0
        f1 = 2 * (precision_metric * sensitivity) / (precision_metric + sensitivity) if (precision_metric + sensitivity) > 0 else 0

        # Append the scores to the lists
        roc_auc_scores.append(roc_auc)
        prc_auc_scores.append(prc_auc)
        sensitivity_scores.append(sensitivity)
        specificity_scores.append(specificity)
        precision_scores.append(precision_metric)
        npv_scores.append(npv)
        f1_scores.append(f1)

    # Calculate confidence intervals
    def calculate_ci(scores):
        lower = np.percentile(scores, (1 - interval) / 2 * 100)
        upper = np.percentile(scores, (1 + interval) / 2 * 100)
        return np.mean(scores), lower, upper

    roc_auc, lower_roc, upper_roc = calculate_ci(roc_auc_scores)
    prc_auc, lower_prc, upper_prc = calculate_ci(prc_auc_scores)
    sensitivity, lower_sens, upper_sens = calculate_ci(sensitivity_scores)
    specificity, lower_spec, upper_spec = calculate_ci(specificity_scores)
    precision_metric, lower_prec, upper_prec = calculate_ci(precision_scores)
    npv, lower_npv, upper_npv = calculate_ci(npv_scores)
    f1, lower_f1, upper_f1 = calculate_ci(f1_scores)

    # Create subplots with gridspec
    fig = plt.figure(figsize=(33/2.54, 10/2.54))  # Figure size in inches, I pass the size in cm and convert to inches
    gs = gridspec.GridSpec(1, 3) 
  
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    gs_right = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[2], height_ratios=[1, 1])  # Ensure equal height for c and d
    ax3_top = fig.add_subplot(gs_right[0])
    ax3_bottom = fig.add_subplot(gs_right[1])

    # Calculate ROC metrics for the original non-bootstrapped results
    fpr, tpr, _ = roc_curve(test_true, test_prob)
    precision, recall, _ = precision_recall_curve(test_true, test_prob)

    # Plot ROC curve
    ax1.plot(fpr, tpr, label=f'AUC = {roc_auc:.3f} \n(95% CI: [{lower_roc:.3f}, {upper_roc:.3f}])')
    ax1.set_xlabel('False Positive Rate', fontsize=14)
    ax1.set_ylabel('True Positive Rate', fontsize=14)
    ax1.set_title('ROC Curve', fontsize=14)
    ax1.legend(loc='lower right')
    ax1.text(-0.1, 1.15, 'a', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

    # Plot PRC curve
    ax2.plot(recall, precision, label=f'AUC = {prc_auc:.3f} \n(95% CI: [{lower_prc:.3f}, {upper_prc:.3f}])')
    ax2.set_xlabel('Recall', fontsize=14)
    ax2.set_ylabel('Precision', fontsize=14)
    ax2.set_title('PRC Curve', fontsize=14)
    ax2.legend(loc='lower left')
    ax2.text(-0.1, 1.15, 'b', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')
    
    # Top table: Confusion matrix of the original non-bootstrapped results
    cm = confusion_matrix(test_true, (test_prob > 0.5).astype(int))
    cm_cell_text = [[f"{cm[i, j]}" for j in range(cm.shape[1])] for i in range(cm.shape[0])]
    cm_cell_text.insert(0, [''] + ['Pred ' + label for label in ['Stroke', 'Seizure']])  # Add column headers
    for i in range(1, len(cm_cell_text)):
        cm_cell_text[i].insert(0, 'True ' + ['Stroke', 'Seizure'][i - 1])  # Add row headers

    ax3_top.axis('off')
    cm_table = ax3_top.table(
        cellText=cm_cell_text,
        loc='center',
        cellLoc='center'
    )
    cm_table.auto_set_font_size(False)
    cm_table.set_fontsize(14)
    cm_table.scale(1.4, 1.6)  # Adjusted height scaling for better text fitting

    # Make row and column labels bold
    for (row, col), cell in cm_table.get_celld().items():
        if row == 0 or col == 0:  # Header cells
            cell.set_text_props(weight='bold')

    ax3_top.text(-0.2, 1.05, 'c', transform=ax3_top.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

    # Bottom table: Classification metrics
    metrics = ['Specificity', 'Sensitivity', 'Precision', 'NPV', 'F1']
    values = [specificity, sensitivity, precision_metric, npv, f1]
    ci_bounds = [(lower_spec, upper_spec), (lower_sens, upper_sens), (lower_prec, upper_prec), (lower_npv, upper_npv), (lower_f1, upper_f1)]
    metrics_cell_text = [[metrics[i], f"{values[i]:.3f}", f"{ci_bounds[i][0]:.3f}-{ci_bounds[i][1]:.3f}"] for i in range(len(metrics))]

    ax3_bottom.axis('off')
    metrics_table = ax3_bottom.table(cellText=metrics_cell_text, colLabels=["Metric", "Score", "95% CI"], loc='center', cellLoc='center')
    metrics_table.auto_set_font_size(False)
    metrics_table.set_fontsize(14)
    metrics_table.scale(1.2, 1.6)

    # Make column headers bold
    for (row, col), cell in metrics_table.get_celld().items():
        if row == 0:  # Header row
            cell.set_text_props(weight='bold')

    ax3_bottom.text(-0.2, 1.05, 'd', transform=ax3_bottom.transAxes, fontsize=14, fontweight='bold', va='top', ha='right')

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.3)
    plt.savefig(f'PMD_extraction_and_analysis/analysis/machine_learning/ml_results.png', dpi=600, bbox_inches='tight')
    plt.show()

    # Print the results    
    print(" ")
    print(f"AUROC = {roc_auc:.4f} (95% CI: [{lower_roc:.4f}, {upper_roc:.4f}])")
    print(f"AUPRC = {prc_auc:.4f} (95% CI: [{lower_prc:.4f}, {upper_prc:.4f}])")
    print(f"Sensitivity = {sensitivity:.4f} (95% CI: [{lower_sens:.4f}, {upper_sens:.4f}])")
    print(f"Specificity = {specificity:.4f} (95% CI: [{lower_spec:.4f}, {upper_spec:.4f}])")
    print(f"Precision = {precision_metric:.4f} (95% CI: [{lower_prec:.4f}, {upper_prec:.4f}])")
    print(f"NPV = {npv:.4f} (95% CI: [{lower_npv:.4f}, {upper_npv:.4f}])")
    print(f"F1 = {f1:.4f} (95% CI: [{lower_f1:.4f}, {upper_f1:.4f}])")
    print(" ")




def metrics_multiclass(test_true, test_prob, nr_bootstraps=5000, interval=0.95):
    """
    Plot AUROC, AUPRC, confusion matrix, and a table with classification metrics side by side in subplots,
    including confidence intervals for AUROC, AUPRC, and other metrics using bootstrapping.

    See https://medium.com/mcd-unison/multiclass-confusion-matrix-clarity-without-confusion-88af1494c1d1
    for a good example of a multiclass confusion matrix. 

    The prc and roc curves are plotted for each class
    """
    class_labels = list(DATASET_PATHS.keys())  # Use class labels from DATASET_PATHS


#--------------------------------------------------------------------------------
# NxN Confusion matrix (Where N is the number of classes) 
#---------------------------------------------------------------------------------
    cm = confusion_matrix(test_true, np.argmax(test_prob, axis=1))
    
    cm_cell_text = [[f"{cm[i, j]}" for j in range(cm.shape[1])] for i in range(cm.shape[0])]
    cm_cell_text.insert(0, [''] + ['Pred ' + label for label in class_labels])  # Add column headers
    for i in range(1, len(cm_cell_text)):
        cm_cell_text[i].insert(0, 'True ' + class_labels[i - 1])  # Add row headers

    # Plot the confusion matrix as a table
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis('off')
    cm_table = ax.table(
        cellText=cm_cell_text,
        loc='center',
        cellLoc='center'
    )
    cm_table.auto_set_font_size(False)
    cm_table.set_fontsize(8)
    cm_table.scale(1.2, 1.2)  # Adjust scaling for better text fitting

    # Make row and column labels bold
    for (row, col), cell in cm_table.get_celld().items():
        if row == 0 or col == 0:  # Header cells
            cell.set_text_props(weight='bold')

    plt.show()


#--------------------------------------------------------------------------------
# Classification metrics per class (with confidence interval) + macro and weighted averages
#---------------------------------------------------------------------------------
    
    # Helper fucntion to calculate the mean and the confidence intervals of the bootstrapped results
    def calculate_ci(scores):
        lower = np.percentile(scores, (1 - interval) / 2 * 100)
        upper = np.percentile(scores, (1 + interval) / 2 * 100)
        return np.mean(scores), lower, upper

    # Initialize a dictionary to store bootstrap results.
    # Each class has a key/value pair, where the key is the class label and the value is a dictionary of metrics.
    # The dictionary of metrics (the sub dictionary). Has a key/value pair for each metric, where the key is the metric name and the value is a list which will be filled with the bootstrap results.
    metrics_bootstrap = {label: {'sensitivity': [], 'specificity': [], 'precision': [], 'npv': [], 'f1': [], 'auroc':[], 'auprc':[]} for label in class_labels}
    
    # Initialize a list to store the mean metrics with confidence intervals
    # macro-averaging calculates the metric per class first, and then averages across the classes.
    # weighted-averaging calculates the metric per class first, and then averages across the classes, but weights each class by the number of samples in that class. Thus classes with more samples have a greater influence on the final metric.
    # micro-averaging calculates the metric globally by counting the total true positives, false negatives, and false positives. Thus classes with more samples have a greater influence on the final metric. But this is not implemented yet.
    metrics_ci = {label: {'sensitivity': [], 'specificity': [], 'precision': [], 'npv': [], 'f1': [], 'auroc':[], 'auprc':[]} for label in class_labels}

    # Bootstrap resampling
    for _ in range(nr_bootstraps):
        indices = np.random.choice(len(test_true), size=len(test_true), replace=True)
        y_true_resampled = test_true[indices]
        y_prob_resampled = test_prob[indices]

        # Calculate metrics for the bootstrapped data
        cm_resampled = confusion_matrix(y_true_resampled, np.argmax(y_prob_resampled, axis=1))
        for i, label in enumerate(class_labels):
            tp = cm_resampled[i, i]
            fn = cm_resampled[i, :].sum() - tp
            fp = cm_resampled[:, i].sum() - tp
            tn = cm_resampled.sum() - (tp + fn + fp)
            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            npv = tn / (tn + fn) if (tn + fn) > 0 else 0
            f1 = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0
            y_true_binary = (y_true_resampled == i).astype(int)
            auroc = roc_auc_score(y_true_binary, y_prob_resampled[:, i])
            auprc = average_precision_score(y_true_binary, y_prob_resampled[:, i])

            # Append the bootstrap results to the metrics_bootstrap dictionary
            metrics_bootstrap[label]['sensitivity'].append(sensitivity)
            metrics_bootstrap[label]['specificity'].append(specificity)
            metrics_bootstrap[label]['precision'].append(precision)
            metrics_bootstrap[label]['npv'].append(npv)
            metrics_bootstrap[label]['f1'].append(f1)
            metrics_bootstrap[label]['auroc'].append(auroc)
            metrics_bootstrap[label]['auprc'].append(auprc)


    # Now that we have all the bootstrap results, we can calculate the mean and the confidence intervals for each metric
    for i, label in enumerate(class_labels):
        sensitivity, lower_sens, upper_sens = calculate_ci(metrics_bootstrap[label]['sensitivity'])
        specificity, lower_spec, upper_spec = calculate_ci(metrics_bootstrap[label]['specificity'])
        precision, lower_prec, upper_prec = calculate_ci(metrics_bootstrap[label]['precision'])
        npv, lower_npv, upper_npv = calculate_ci(metrics_bootstrap[label]['npv'])
        f1, lower_f1, upper_f1 = calculate_ci(metrics_bootstrap[label]['f1'])
        auroc, lower_roc, upper_roc = calculate_ci(metrics_bootstrap[label]['auroc'])
        auprc, lower_prc, upper_prc = calculate_ci(metrics_bootstrap[label]['auprc'])

        # Add the mean + ci of the metrics to the metrics_ci dictionary
        metrics_ci[label] = {
            'sensitivity': [sensitivity, lower_sens, upper_sens],
            'specificity': [specificity, lower_spec, upper_spec],
            'precision': [precision, lower_prec, upper_prec],
            'npv': [npv, lower_npv, upper_npv],
            'f1': [f1, lower_f1, upper_f1],
            'auroc': [auroc, lower_roc, upper_roc],
            'auprc': [auprc, lower_prc, upper_prc]
        }

    # Now that we have the mean and the confidence intervals for each metric and class, we can calculate the macro and weighted averages
    # Calculate the macro average of the metrics 
    macro_avg_sensitivity = np.mean([metrics_ci[label]['sensitivity'][0] for label in class_labels])
    macro_avg_specificity = np.mean([metrics_ci[label]['specificity'][0] for label in class_labels])
    macro_avg_precision = np.mean([metrics_ci[label]['precision'][0] for label in class_labels])   
    macro_avg_npv = np.mean([metrics_ci[label]['npv'][0] for label in class_labels])
    macro_avg_f1 = np.mean([metrics_ci[label]['f1'][0] for label in class_labels])
    macro_avg_auroc = np.mean([metrics_ci[label]['auroc'][0] for label in class_labels])
    macro_avg_auprc = np.mean([metrics_ci[label]['auprc'][0] for label in class_labels])
    
    # Calculate the weighted average of the metrics
    class_weights = {label: (test_true == i).sum() / len(test_true) for i, label in enumerate(class_labels)}
    weighted_avg_sensitivity = np.sum([metrics_ci[label]['sensitivity'][0] * class_weights[label] for label in class_labels])
    weighted_avg_specificity = np.sum([metrics_ci[label]['specificity'][0] * class_weights[label] for label in class_labels])
    weighted_avg_precision = np.sum([metrics_ci[label]['precision'][0] * class_weights[label] for label in class_labels])
    weighted_avg_npv = np.sum([metrics_ci[label]['npv'][0] * class_weights[label] for label in class_labels])
    weighted_avg_f1 = np.sum([metrics_ci[label]['f1'][0] * class_weights[label] for label in class_labels])
    weighted_avg_auroc = np.sum([metrics_ci[label]['auroc'][0] * class_weights[label] for label in class_labels])
    weighted_avg_auprc = np.sum([metrics_ci[label]['auprc'][0] * class_weights[label] for label in class_labels])


    # Add the macro and weighted averages to the metrics_ci dictionary
    metrics_ci['macro avg'] = {
        'sensitivity': [macro_avg_sensitivity],
        'specificity': [macro_avg_specificity],
        'precision': [macro_avg_precision],
        'npv': [macro_avg_npv],
        'f1': [macro_avg_f1],
        'auroc': [macro_avg_auroc],
        'auprc': [macro_avg_auprc]
    }	

    metrics_ci['weighted avg'] = {
        'sensitivity': [weighted_avg_sensitivity],
        'specificity': [weighted_avg_specificity],
        'precision': [weighted_avg_precision],
        'npv': [weighted_avg_npv],
        'f1': [weighted_avg_f1],
        'auroc': [weighted_avg_auroc],
        'auprc': [weighted_avg_auprc]
    }
        
    # Prepare data for the table
    metrics_table_data = []
    for label, metrics in metrics_ci.items():
        row = [label]
        for metric in ['sensitivity', 'specificity', 'precision', 'npv', 'f1', 'auroc', 'auprc']:
            if len(metrics[metric]) == 3:  # Mean, lower CI, upper CI
                mean, lower, upper = metrics[metric]
                row.append(f"{mean:.3f} ({lower:.3f}-{upper:.3f})")
            else:  # For macro and weighted averages without CI
                row.append(f"{metrics[metric][0]:.3f}")
        metrics_table_data.append(row)

    # Create a second plot for classification metrics per class
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.axis('off')
    metrics_table = ax.table(
        cellText=metrics_table_data,
        colLabels=["Class", "Sensitivity (95% CI)", "Specificity (95% CI)", "Precision (95% CI)", "NPV (95% CI)", "F1 (95% CI)", "AUROC (95% CI)", "AUPRC (95% CI)"],
        loc='center',
        cellLoc='center'
    )

    metrics_table.auto_set_font_size(False)
    metrics_table.set_fontsize(8)
    metrics_table.scale(1.2, 1.2)  # Adjust scaling for better text fitting

    # Make column headers bold
    for (row, col), cell in metrics_table.get_celld().items():
        if row == 0:  # Header row
            cell.set_text_props(weight='bold')

    plt.title("Classification Metrics Per Class with Confidence Intervals", fontsize=14, fontweight='bold')
    plt.show()



#--------------------------------------------------------------------------------
# ROC and PRC curves for each class
#---------------------------------------------------------------------------------

    for i, class_label in enumerate(class_labels):
        # Binarize the true labels for the current class
        test_true_binary = (test_true == i).astype(int)

        # Calculate ROC metrics
        fpr, tpr, _ = roc_curve(test_true_binary, test_prob[:, i])
        roc_auc = roc_auc_score(test_true_binary, test_prob[:, i])

        # Calculate PRC metrics
        precision, recall, _ = precision_recall_curve(test_true_binary, test_prob[:, i])
        prc_auc = average_precision_score(test_true_binary, test_prob[:, i])

        # Create a new figure with 2 subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

        # Plot ROC curve
        ax1.plot(fpr, tpr, label=f'AUC = {roc_auc:.3f}')
        ax1.set_xlabel('False Positive Rate', fontsize=12)
        ax1.set_ylabel('True Positive Rate', fontsize=12)
        ax1.set_title(f'ROC Curve - {class_label}', fontsize=14)
        ax1.legend(loc='lower right')

        # Plot PRC curve
        ax2.plot(recall, precision, label=f'AUC = {prc_auc:.3f}')
        ax2.set_xlabel('Recall', fontsize=12)
        ax2.set_ylabel('Precision', fontsize=12)
        ax2.set_title(f'PRC Curve - {class_label}', fontsize=14)
        ax2.legend(loc='lower left')

        # Adjust layout and show the figure
        plt.tight_layout()
        plt.show()