import argparse
import textwrap
import pandas as pd
import numpy as np
from scipy import stats


def _add_x_input_to_parser(parser):
    parser.add_argument(
        '--x',
        dest='X',
        help='path to input X dataframe (pickle only)')

    
def _add_y_input_to_parser(parser):
    parser.add_argument(
        '--y',
        dest='y',
        help='path to input y dataframe (pickle only)')


def _check_for_missing_values(df):
    return df.isnull().values.any()


def t_test(t_val, df):
    return (1 - stats.t.cdf(abs(t_val), df=df)) * 2


def _create_parser():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent('''\
            Hypothesis testing using linear regression theory. '''),
        epilog=textwrap.dedent('''\
        Examples:             
        python linear_model_hypothesis_testing.py
        --x _sample_data/x.pkl --y _sample_data/y.pkl ''')
    )

    _add_x_input_to_parser(parser)
    _add_y_input_to_parser(parser)

    return parser


def _hypothesis_testing(X, y):
    XtX = X.transpose().dot(X)
    XtX_inv = pd.DataFrame(np.linalg.pinv(XtX.values), XtX.columns, XtX.index)
    
    beta_hat = XtX_inv.dot(X.transpose()).dot(y)

    n = X.shape[0]  # num of samples
    k = X.shape[1]  # num of columns (features)
    p = k + 1

    SSres = y.dot(y.transpose()) - beta_hat.transpose().dot(X.transpose()).dot(y)
    sigma_hat = np.sqrt(SSres / (n - p))

    C = XtX_inv
    sqrtCii = pd.DataFrame(np.sqrt(np.diag(C)), index=[C.columns], columns=['sqrtCii'])

    t = (beta_hat - 0).divide(sigma_hat * sqrtCii['sqrtCii'])
    p_values = t.apply(lambda t_val: t_test(t_val, df=(n - p)))
    
    print '\nP-VALUES RESULTS:\n------------------\n'
    print p_values

    
def main():
    parser = _create_parser()
    args = parser.parse_args()

    X = pd.read_pickle(args.X)
    y = pd.read_pickle(args.y)
    
    if _check_for_missing_values(X) or _check_for_missing_values(y):
        raise ValueError('Inputs cannot contain NaN (missing values))')

    _hypothesis_testing(X, y)

    
if __name__ == '__main__':
    main()
