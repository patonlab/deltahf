"""Tests for CLI entry point."""

import pytest

from deltahf.__main__ import build_parser, main


class TestParser:
    def test_help_exits(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            build_parser().parse_args(["--help"])
        assert exc_info.value.code == 0

    def test_fit_subcommand(self):
        args = build_parser().parse_args(["fit", "-i", "training.csv"])
        assert args.command == "fit"
        assert args.input == "training.csv"

    def test_predict_subcommand(self):
        args = build_parser().parse_args(["predict", "-i", "molecules.csv", "--epsilon", "fitted.json"])
        assert args.command == "predict"
        assert args.input == "molecules.csv"
        assert args.epsilon == "fitted.json"

    def test_fit_defaults(self):
        args = build_parser().parse_args(["fit", "-i", "training.csv"])
        assert args.model == "both"
        assert args.kfold == 10
        assert args.n_conformers == 5

    def test_predict_defaults(self):
        args = build_parser().parse_args(["predict", "-i", "mol.csv", "--epsilon", "e.json"])
        assert args.model == "4param"
        assert args.n_conformers == 5

    def test_no_subcommand_prints_help(self, capsys):
        with pytest.raises(SystemExit):
            main([])
