from prep_data import prepare_tsfresh_train_test
import os

if __name__ == "__main__":
    ROOT = os.path.dirname(os.path.abspath(__file__))
    tsfresh_path = "feature_small.csv"
    labels_path = "label/labels_small.csv"

    X_train_sel, X_test_sel, y_train, y_test = prepare_tsfresh_train_test(ROOT,
    tsfresh_path,
    labels_path,
    target_col="Maincluster",
    test_size=0.2,
    random_state=42,
    use_relevant_features=True,
    fdr_level=0.05,
    n_jobs=0,
)
