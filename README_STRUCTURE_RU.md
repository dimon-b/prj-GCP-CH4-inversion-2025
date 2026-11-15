# GCP/GMB 2025 — обновлённые файлы и структура

Распаковать **содержимое архива** в существующий `prj_GCP_v25/`, чтобы получился каталог `prj_GCP_v25/c_gcpv3/`. При обновлении проекта перенести только показанные здесь файлы; чужие/старые файлы проекта не удаляются.

```
prj_GCP_v25/
├── c_gcpv3/
│   ├── _main.py
│   ├── _set_case.py
│   ├── test/
│   │   ├── test_all.py, test_config.py, test_f2.py ...
│   │   └── reference/legacy_main_actmstn.f90   (только эталон для теста F6)
│   ├── fortran/
│   │   ├── f0_*.f90, f2_*.f90, f3_*.f90 ... f7_*.f90
│   │   ├── lsfbf4_main.f, lsfbf4_sub.f         (неизменённая библиотека Nakazawa)
│   │   ├── build_fortran.sh
│   │   └── bin/                                 (создаётся после компиляции)
│   └── inversion/
│       ├── __init__.py
│       ├── inversion_runner.py
│       ├── check_inv_inputs.py
│       └── control.example
├── inv_dir/
├── fluxes/
├── trout_MIROC/
└── results2025/
```

## Что изменено в `_set_case.py`

- `self.fort_dir = './fortran/'`
- `self.fort_exe_dir = './fortran/bin/'` по умолчанию (доступен `CH4_FORTRAN_BIN_DIR`)
- `self.inv_code_dir = './inversion/'`
- `self.inv_control_template = './inversion/control.example'` по умолчанию
- `self.exe_inv = './inversion/tdi_ch4'` по умолчанию (можно переопределить через `CH4_EXE_INV`)

**Все остальные каталоги и пользовательские пути из объединённого `_set_case.py` оставлены без изменений.** Команды запускать из `prj_GCP_v25/c_gcpv3/`. Для всех относительных путей `../inv_dir/`, `../results2025/`, `../trout_MIROC/` это обязательно.

Исходный executable инверсии `tdi_ch4` в архив **не добавлялся**: мы переносим только ранее отредактированные исходники. Для реального запуска инверсии надо скопировать или сделать симлинк из вашего существующего `codeint54/tdi_ch4` в `c_gcpv3/inversion/tdi_ch4` либо выставить `CH4_EXE_INV`. Исходные Fortran-программы solver и вспомогательные файлы `codeint54/files4run/` не изменены и не перемещались. В `control.example` ссылки на них остались исходными: проверьте доступность на сервере.

## Компиляция и тесты

```bash
cd /S/data01/G5070/y0715/prj_GCP_v25/c_gcpv3
bash fortran/build_fortran.sh
python test/test_all.py
python test/test_f3_binary.py
python _main.py --case CYC --steps f2 f3 f6 f7 --dry-run
python _main.py --case INCA --steps f2 f3 f6 f7 --dry-run
python _main.py --case CYC --steps inversion --dry-run
```

`test/test_f3_binary.py` выполняется отдельно после компиляции. Никакие реальные NetCDF или GT3 не используются тестами без явного запуска основных этапов.

Версии F0–F7 и численный алгоритм фильтра Накадзавы **не изменялись** при этой реорганизации; исправлены только импорты, ссылки тестов на файлы и каталоги executable в `_set_case.py`.
