# Литературная проверка новизны

Вопрос: публиковал ли кто-нибудь наблюдаемый скейлинг
`λ₂(normalized Laplacian, giant) ~ n^(−α)` при α ≈ 0.15…0.17
для power-law конфигурационной модели с **min degree 2** (d_max = 50).

Метод: arXiv API (полнотекстовые запросы) + Crossref. Semantic Scholar отдавал 429.
Дата: 2026-09-23.

---

## Итог (предварительно)

**Прямого совпадения не найдено.** Никто, судя по доступному поиску, не публиковал
именно этот скейлинг. Это **предварительный** вывод: полнотекстовый поиск arXiv
слабый, citation-граф и MathSciNet/Web of Science не проверялись.

## Ключевой аргумент в пользу новизны

Все теоремы, дающие gap = Θ(1), требуют одного из двух:

- **min degree ≥ 3** (Hermon–Li–Yao–Zhang), либо
- **d̄_min ≫ ln²n** (Chung–Lu–Vu).

Наш случай (heavy tail + **min degree 2**, d̄_min = O(1)) лежит **ровно в
запрещённой зоне** — там, где теория не даёт Θ(1). Наблюдаемый `n^(−0.15)`
заполняет именно этот пробел.

## Ближайшие работы

| Работа | Что доказывает | Отношение к нам |
|---|---|---|
| **Hermon, Li, Yao, Zhang**, arXiv:2105.11585, *Ann. Prob.* 50(5) 2022 | config model со степенями в **[3, d̄]**, spectral gap / relaxation time | **ключевое**: требуют min degree ≥3 — ровно наш «const»-режим; d_min=2 вне условий |
| **Chung–Lu–Vu** (semicircle) | концентрация при **d̄_min ≫ ln²n** | наш режим d̄_min=O(1) — вне |
| **Wang, Li, Huang**, arXiv:2412.02087 (2024) | semicircle для норм. лапласиана config model, «mild assumptions» | про bulk-спектр, не про gap; heavy-tail+deg2, похоже, не покрыт |
| **Coste**, arXiv:1708.00530 (2017) | spectral gap разреженных random **digraph** по degree sequence (trace method) | методология |
| **Cai, Caputo, Perarnau, Quattropani**, arXiv:2104.08389 | directed CM, heavy-tailed in-degrees, mixing, power-law stationary | контекст heavy-tail |
| **Fernley**, arXiv:2211.13537 | voter model на scale-free, mixing, τ∈(2,3] | контекст «ultrasmall world» |
| **Lu & Peng**, arXiv:1204.6207 | Laplacian concentration требует δ ≫ ln n | почему sparse тяжёл |
| **Ducatez & Rivier**, EJP 2025, 10.1214/25-ejp1303 | gap лапласиана ER и вложенные деревья | ER-giant — экспандер; контраст с нами |
| **Avrachenkov et al.**, arXiv:1910.08869 | smallest nonzero eigenvalue → 0, power-law плотность у 0, spectral dimension | похожая феноменология |

## Что НЕ покрыто поиском (риски)

- Citation-граф: кто цитирует Coja-Oghlan–Lanka, Hermon et al., Bordenave — надёжнее всего, но arXiv API этого не даёт.
- Google Scholar / Web of Science / MathSciNet.
- Запросы про «spectral dimension» scale-free могли пропустить уже известный показатель.

## Что нужно для надёжной проверки

1. Targeted-поиск «кто цитирует» по 2–3 ключевым работам (Crossref citations / Semantic Scholar с паузой).
2. Проверить Google Scholar / MathSciNet.
3. Тест «универсальности» α: варьировать `d_max` (фиксирован d_max=50; если α сдвинется — эффект обрезки, а не хвоста).
