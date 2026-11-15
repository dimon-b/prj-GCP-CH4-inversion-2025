# INCA / CYC — единое обновление подготовки инверсии

**Применять только после завершения текущего расчёта.**
Текущую папку `results2025/...` не менять и не удалять.

Это один комплект, объединяющий два обновления:

1. **Input snapshots**: вместо символических ссылок — независимые копии входных файлов в каждом новом расчёте.
2. **No vardump**: отключена запись большой полной матрицы posterior covariance, при сохранении расчёта SVD и posterior fluxes.

## Какие файлы заменить

Из корня проекта `/S/data01/G5070/y0715/prj_GCP_v25/` распакуйте архив с сохранением путей:

```text
c_gcpv3/inversion/inversion_runner.py
c_gcpv3/inversion/control.example
c_gcpv3/test/test_inversion.py
c_gcpv3/inversion/README_UPDATE_RU.md
```

Если есть собственные изменения в `control.example`, сравните их до замены. Поставляемый шаблон соответствует текущей конфигурации: `constf='files4run/spatial-consts-none'`, период 1999–2025, путь `obsrvCH4_60sites.txt` (57 активных станций в его содержимом), `dumpmats=.false.`.

Не нужно менять `_main.py`, `_set_case.py`, `tdi.f`, `read_globalview.f90`, `Makefile` или перекомпилировать `tdi_ch4`.

## Как подготавливаются новые расчёты

Для каждого **нового** эксперимента в `results2025/CASE/icase_prior/icase_unc_model/` создаётся:

```text
control                   # использует локальный obsrvCH4_sites.txt
obsrvCH4_9925.nc          # реальная копия F7, не symlink
sourcepriortdi.nc         # реальная копия F5, не symlink
obsrvCH4_sites.txt        # реальная копия списка станций
files4run/                # реальная копия служебных файлов
  spatial-consts-none
  spatial-srccor
  spatial-groups-tdi
  new-recipes
  ...
tdi_ch4                   # сохранённый executable этой инверсии
RUN_INPUTS.txt            # источники, SHA-256, параметры, сводка Green's functions
RUN_INPUTS.sha256         # контрольные суммы локальных копий
prior.sources
predicted.sources
...
```

`inversion_runner.py` вычисляет SHA-256 до/после копирования и проверяет копию. Если источник поменялся во время подготовки, запуск прекращается с `SNAPSHOT_FAILED.txt`.

Каталог результата **не перезаписывается**, если уже существует. Для нового эксперимента нужен новый `icase` или отдельное заранее подготовленное имя (не удаляйте текущий каталог расчёта).

### Green's functions

`greenf` оставлен ссылкой на общий абсолютный **путь** `inv_dir/gmats/` в control (из-за объёма). Этот каталог **не копируется**. В `RUN_INPUTS.txt` сохраняются путь, число файлов, суммарный размер и SHA-256 инвентаря (имена, размеры, mtime), но **не хеш содержимого GF**. Не меняйте `gmats/` во время расчёта; после его завершения скрипт сравнивает инвентарь повторно.

### Отключён `vardump`

В `control.example`:

```fortran
dumpmats = .false.
```

Даже если старый шаблон содержит `dumpmats=.true.`, runner принудительно запишет `.false.` в **новый локальный** `control`. Файл `vardump` не создаётся. Решение Bayesian inversion, SVD и `prior.sources`/`predicted.sources` от этого не меняются.

Старые `vardump` не удаляются автоматически. `RUN_INPUTS.txt` содержит `vardump_enabled=false`.

## Проверка после замены

Из `c_gcpv3/`:

```bash
python test/test_inversion.py
```

## Запуск нового эксперимента

После завершения текущей задачи и выбора **нового `icase`**, например:

```bash
nohup python -u _main.py --case INCA --steps inversion --inv-prior p30 --inv-unc ctl > inversion_INCA_next.log 2>&1 < /dev/null &
```

Проверка локальных файлов нового расчёта:

```bash
cd ../results2025/INCA/<icase>_p30/<icase>_ctl_inv1
ls -l obsrvCH4_9925.nc sourcepriortdi.nc obsrvCH4_sites.txt files4run tdi_ch4
sha256sum -c RUN_INPUTS.sha256
grep -i '^dumpmats' control
ls vardump
```

Ожидаются: обычные файлы (без `->`), `OK` для всех SHA-256, `dumpmats = .false.`, а `ls vardump` сообщает об отсутствии файла.

`self.icase` пока задаёт имя эксперимента, **но не выбирает сеть наблюдений**. Автоматическое соответствие `s057 → obsrvCH4_57sites.txt` остаётся в TODO.
