#!/usr/bin/env python

import os
import pandas as pd

def read_sacct(path):

    df = pd.read_csv(path, sep='|')

    return df

def clean_up(df, cols):

    df = df.drop(cols, 1, inplace=True)

    return df

def simplify_ix(df):
    df = df[(df['JobName'] != 'batch') & (df['JobName'] != 'python')]
    return df

def print_full(df):
    pd.set_option('display.madf_rows', len(df))
    print(df)
    pd.reset_option('display.madf_rows')

def incomplete(df):

    df = df[df.State != 'COMPLETED']

    return df


