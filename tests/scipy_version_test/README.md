# scipy Version Test

Tests whether `scipy.signal.savgol_filter` behavior changed between scipy 1.4.1 (Jan 2020) and scipy 1.16.3 (2025).

```bash
pixi run --environment scipy14 python test_savgol.py > scipy14_output.txt
pixi run --environment scipy116 python test_savgol.py > scipy116_output.txt
diff scipy14_output.txt scipy116_output.txt
```

**Difference: 1.4×10⁻¹³** (14th decimal place - floating point noise)
